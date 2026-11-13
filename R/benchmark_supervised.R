#!/usr/bin/env Rscript
# Dos metodos publicados de prediccion de hora circadiana, ZeitZeiger (Hughey 2016, 2017) y tauFisher
# (Duan, Ngo et al. 2024), sobre los mismos datos y los mismos pliegues que el estimador del articulo.
# Lo llama python/analyses/32_benchmark_methods.py, que escribe las entradas y lee las salidas.
#
#   Rscript R/benchmark_supervised.R <directorio de trabajo>
#
# Entradas, en el directorio de trabajo:
#   train_expr.tsv   genes x muestras, log2 CPM, cohorte de entrenamiento con horas registradas
#   train_meta.tsv   sample, subject, clock_h
#   test_expr.tsv    genes x muestras de la cohorte externa (otro tejido, otra plataforma)
#   test_meta.tsv    sample, subject, clock_h
#   clock_genes.txt  genes del reloj, del registro de genes del repositorio
# Salidas: predictions_loso.tsv y predictions_transfer.tsv, con la hora predicha por cada metodo.
#
# Cada metodo se usa como lo documentan sus autores, sin ajustar nada con las muestras de prueba:
#   ZeitZeiger  sumabsv = 2, nSpc = 2, los valores de sus articulos
#   tauFisher   genes periodicos por JTK (MetaCycle) en el entrenamiento, los 20 mas significativos
#               mas los genes del reloj; curvas FDA con 5 funciones base; regresion multinomial
suppressMessages({library(zeitzeiger); library(tauFisher); library(MetaCycle)})
args <- commandArgs(trailingOnly = TRUE); wd <- args[1]
rd <- function(f) { x <- read.delim(file.path(wd, f), check.names = FALSE, stringsAsFactors = FALSE); x }
E  <- rd("train_expr.tsv"); rownames(E) <- E[[1]]; E <- as.matrix(E[, -1])
M  <- rd("train_meta.tsv"); stopifnot(all(M$sample == colnames(E)))
E2 <- rd("test_expr.tsv");  rownames(E2) <- E2[[1]]; E2 <- as.matrix(E2[, -1])
M2 <- rd("test_meta.tsv");  stopifnot(all(M2$sample == colnames(E2)))
CLOCK <- toupper(readLines(file.path(wd, "clock_genes.txt")))      # del registro de genes del repositorio

zz <- function(xtr, ttr, xte) {
  fit <- zeitzeiger(t(xtr), (ttr %% 24) / 24, t(xte), sumabsv = 2, nSpc = 2)
  as.numeric(fit$predResult$timePred) * 24
}

periodic <- function(xtr, ttr) {
  # JTK sobre la serie de entrenamiento (media por hora), como en la vineta de tauFisher
  tt <- sort(unique(ttr)); avg <- sapply(tt, function(h) rowMeans(xtr[, ttr == h, drop = FALSE]))
  d <- data.frame(Gene = rownames(xtr), avg, check.names = FALSE); colnames(d)[-1] <- tt
  f <- tempfile(fileext = ".txt"); o <- tempfile(); dir.create(o)
  write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(capture.output(suppressWarnings(meta2d(infile = f, filestyle = "txt", outdir = o, timepoints = tt,
      cycMethod = "JTK", minper = 24, maxper = 24, outIntegration = "noIntegration", outputFile = TRUE))))
  r <- read.delim(list.files(o, pattern = "^JTKresult", full.names = TRUE)[1])
  r <- r[order(r$ADJ.P, -r$AMP), ]
  unique(toupper(c(head(r$CycID, 20), intersect(CLOCK, rownames(xtr)))))
}

tf <- function(xtr, ttr, xte) {
  out <- tryCatch({
    g <- periodic(xtr, ttr)
    o <- order(ttr); lin_tr <- 2^xtr[, o] - 1; lin_tr[lin_tr < 0] <- 0; rownames(lin_tr) <- toupper(rownames(xtr))
    nrep <- max(table(ttr))
    invisible(capture.output(m <- train_tauFisher(g, lin_tr, ttr[o], nrep = nrep, numbasis = 5)))
    lin_te <- 2^xte - 1; lin_te[lin_te < 0] <- 0; rownames(lin_te) <- toupper(rownames(xte))
    keep <- m$Genes[m$Genes %in% rownames(lin_te)]
    if (length(keep) < length(m$Genes)) stop("faltan genes del modelo en la cohorte de prueba")
    invisible(capture.output(p <- test_tauFisher(m$Genes, lin_te, rep(0, ncol(lin_te)), m$PCA_embedding, m$Model)))
    as.numeric(as.character(p$predicted))
  }, error = function(e) { message("tauFisher: ", conditionMessage(e)); rep(NA_real_, ncol(xte)) })
  out
}

# ---- 1. dejando sujetos enteros fuera, en la cohorte de entrenamiento
res <- list()
for (s in unique(M$subject)) {
  te <- M$subject == s
  res[[length(res) + 1]] <- data.frame(sample = M$sample[te], subject = s, clock_h = M$clock_h[te],
      zeitzeiger_h = zz(E[, !te], M$clock_h[!te], E[, te, drop = FALSE]),
      taufisher_h  = tf(E[, !te], M$clock_h[!te], E[, te, drop = FALSE]))
  cat("  sujeto", s, "listo\n")
}
write.table(do.call(rbind, res), file.path(wd, "predictions_loso.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# ---- 2. transferencia: entrenar en una cohorte y predecir otra, de otro tejido y otra plataforma
common <- intersect(rownames(E), rownames(E2))
# Entre plataformas la escala de cada gen cambia. ZeitZeiger se prueba ademas con cada gen
# estandarizado dentro de su cohorte, que no usa ninguna hora de la cohorte de prueba. (Sus autores
# recomiendan ComBat para esto; exige el paquete sva y no cambia que la correccion sea sin etiquetas.)
zs <- function(x) { s <- apply(x, 1, sd); x <- x[s > 0, , drop = FALSE]; t(scale(t(x))) }
A <- zs(E[common, ]); B <- zs(E2[common, ]); cz <- intersect(rownames(A), rownames(B))
tr <- data.frame(sample = M2$sample, subject = M2$subject, clock_h = M2$clock_h,
    zeitzeiger_h = zz(E[common, ], M$clock_h, E2[common, ]),
    zeitzeiger_scaled_h = zz(A[cz, ], M$clock_h, B[cz, ]),
    taufisher_h  = tf(E[common, ], M$clock_h, E2[common, ]))
write.table(tr, file.path(wd, "predictions_transfer.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("ok\n")
