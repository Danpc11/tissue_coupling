#!/usr/bin/env bash
# Descarga los datos publicos que espera tissue_coupling (docs/DATA.md, 00_prepare_data.py)
# en data/raw: FTP de GEO, GCS de GTEx, Broad (MSigDB), OmniPath, recount3, ARCHS4 (opcional) y
# www.ncbi.nlm.nih.gov para las cuentas generadas por NCBI. Este ultimo puede responder con una
# pagina de reCAPTCHA a clientes no interactivos; en ese caso se avisa y hay que reintentar (las
# alternativas de respaldo estan documentadas en docs/DATA.md).
#
#   bash download_public.sh            # todo menos ARCHS4
#   T2D_ARCHS4_DL=1 bash download_public.sh   # ademas el HDF5 de ARCHS4 (62 GB) para GSE159984
set -uo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
RAW=$ROOT/data/raw
UA='Mozilla/5.0'
mkdir -p "$RAW/msigdb" "$RAW/_soft"
cd "$RAW"
get () {  # get <url> <destino>
  local url=$1 dst=$2
  if [ -s "$dst" ]; then echo "  [ya]    $dst"; return 0; fi
  if curl -s -S -f -A "$UA" --retry 5 --retry-delay 10 -L -o "$dst.part" "$url"; then
    mv "$dst.part" "$dst"; echo "  [ok]    $dst ($(du -h "$dst" | cut -f1))"
  else
    rm -f "$dst.part"; echo "  [FALLO] $dst  <- $url"; return 1
  fi
}
listdir () { curl -s -A "$UA" "$1" | grep -o 'href="[^"]*"' | sed 's/href="//; s/"$//' | grep -v '^\.\./\|^/\|^http'; }

echo "### $(date)  inicio"
echo "## 1. Matrices de serie (GEO FTP)"
for a in GSE13070 GSE17992 GSE25462 GSE26637 GSE41168 GSE66446 GSE78721 GSE87761 \
         GSE40589 GSE95674 GSE95675 GSE107537 GSE108539 GSE129843 GSE148640 GSE159984 GSE182117 GSE182120; do
  n=${a:0:$((${#a}-3))}nnn
  files=$(listdir "https://ftp.ncbi.nlm.nih.gov/geo/series/$n/$a/matrix/" | grep '_series_matrix.txt.gz$')
  [ -z "$files" ] && echo "  [FALLO] $a: sin matrices listadas"
  for f in $files; do get "https://ftp.ncbi.nlm.nih.gov/geo/series/$n/$a/matrix/$f" "$f"; done
  sleep 1
done

echo "## 2. Suplementarios de autor (cuentas y procesados)"
while read -r a f; do
  n=${a:0:$((${#a}-3))}nnn
  get "https://ftp.ncbi.nlm.nih.gov/geo/series/$n/$a/suppl/$f" "$f"
done <<'LIST'
GSE108539 GSE108539_norm_reads_count.txt.gz
GSE148640 GSE148640_processed.txt.gz
GSE182117 GSE182117_counts.tsv.gz
GSE182117 GSE182117_logCPM_batchCorrected.tsv.gz
GSE129843 GSE129843_README.txt
GSE129843 GSE129843_RESTRICT.txt.gz
GSE148640 filelist.txt
LIST
mv -f filelist.txt GSE148640_filelist.txt 2>/dev/null
for f in GSE159984_de-t2d.tsv.gz; do get "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE159nnn/GSE159984/suppl/$f" "$f"; done

echo "## 3. Anotacion de plataformas (tabla del SOFT, leida en streaming hasta !platform_table_end)"
# Los *_family.soft.gz incluyen todas las muestras de la plataforma (GPL570: decenas de GB,
# GPL15207: 4.6 GB), pero la tabla de la plataforma va al principio: se descarga solo hasta ahi.
soft_table () {  # soft_table <url> <destino>
  [ -s "$2" ] && { echo "  [ya]    $2"; return 0; }
  curl -s -A "$UA" "$1" 2>/dev/null | zcat 2>/dev/null \
    | awk '/^!platform_table_begin/{f=1;next} /^!platform_table_end/{exit} f' > "$2.part"
  if [ -s "$2.part" ] && head -1 "$2.part" | grep -q '^ID'; then
    mv "$2.part" "$2"; echo "  [ok]    $2 ($(wc -l < "$2") lineas)"
  else rm -f "$2.part"; echo "  [FALLO] $2 <- $1"; fi
}
soft_table https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/soft/GPL570_family.soft.gz     GPL570_annot.txt
soft_table https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL15nnn/GPL15207/soft/GPL15207_family.soft.gz GPL15207_annot.txt
soft_table https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL24nnn/GPL24324/soft/GPL24324_family.soft.gz GPL24324-69524.txt
get "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz" "_soft/GPL570.annot.gz"

echo "## 4. GTEx v11 (GCS publico)"
B=https://storage.googleapis.com/adult-gtex
get $B/annotations/v11/metadata-files/GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt GTEx_Analysis_v11_Annotations_SampleAttributesDS.txt
get $B/annotations/v11/metadata-files/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt GTEx_Analysis_v11_Annotations_SubjectPhenotypesDS.txt
get $B/annotations/v11/metadata-files/GTEx_Analysis_v11_Annotations_SampleAttributesDD.xlsx GTEx_Analysis_v11_Annotations_SampleAttributesDD.xlsx
get $B/annotations/v11/metadata-files/GTEx_Analysis_v11_Annotations_SubjectPhenotypesDD.xlsx GTEx_Analysis_v11_Annotations_SubjectPhenotypesDD.xlsx
for t in muscle_skeletal adipose_subcutaneous adipose_visceral_omentum whole_blood; do
  get $B/bulk-gex/v11/rna-seq/counts-by-tissue/gene_reads_v11_$t.gct.gz gene_reads_v11_$t.gct.gz
done

echo "## 5. MSigDB 2026.2 (H, C2, C5, C5-HPO, C8; simbolos)"
M=https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.2.Hs
for c in h.all c2.all c5.all c5.hpo c8.all; do
  get $M/$c.v2026.2.Hs.symbols.gmt msigdb/$c.v2026.2.Hs.symbols.gmt
done

echo "## 6. DoRothEA A-C (OmniPath)"
get 'https://omnipathdb.org/interactions?datasets=dorothea&dorothea_levels=A,B,C&organisms=9606&genesymbols=yes&fields=dorothea_level,sources,references' dorothea_ABC.tsv

echo "## 7. Anotacion de genes NCBI (sustituto de Human.GRCh38.p13.annot.tsv.gz)"
get https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz _soft/Homo_sapiens.gene_info.gz
if [ ! -s Human_GRCh38_p13_annot_tsv.gz ]; then
  # columnas que leen los scripts: GeneID y Symbol (gene_info las trae con ese mismo nombre)
  zcat _soft/Homo_sapiens.gene_info.gz | awk -F'\t' 'BEGIN{OFS="\t"; print "GeneID","Symbol","Description","Synonyms","GeneType"} NR>1{print $2,$3,$9,$5,$10}' | gzip > Human_GRCh38_p13_annot_tsv.gz.part && mv Human_GRCh38_p13_annot_tsv.gz.part Human_GRCh38_p13_annot_tsv.gz
  echo "  [ok]    Human_GRCh38_p13_annot_tsv.gz (DERIVADO de gene_info; ver README_descarga.md)"
fi
ln -sfn Human_GRCh38_p13_annot_tsv.gz Human_GRCh38_p13_annot.tsv.gz

echo "## 8. Otros ficheros que leen 00b_build_matrices.py y 00c_archs4_counts.py"
soft_table https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL17nnn/GPL17586/soft/GPL17586_family.soft.gz GPL17586_annot.txt
soft_table https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL4nnn/GPL4133/soft/GPL4133_family.soft.gz     GPL4133_annot.txt
# GPL6848 (Agilent 4x44K, nombres de sonda): GSE40589, clamp frente a salina a las mismas horas.
# No tiene fichero .annot curado en GEO; la tabla se toma del SOFT como las demas.
soft_table https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL6nnn/GPL6848/soft/GPL6848_family.soft.gz     GPL6848_annot.txt
soft_table https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL15nnn/GPL15331/soft/GPL15331_family.soft.gz GPL15331_annot.txt
for c in c2.cp.reactome c2.cp.wikipathways c2.cp.kegg_medicus c2.cgp c2.cp.pid c3.tft.gtrd c3.mir.mirdb c4.all c5.go.bp c7.all; do
  get $M/$c.v2026.2.Hs.symbols.gmt msigdb/$c.v2026.2.Hs.symbols.gmt
done
# cuentas generadas por NCBI (www.ncbi.nlm.nih.gov; si responde con una pagina HTML, reintentar mas tarde)
N='https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts'
# El CAPTCHA de NCBI es intermitente: se reintenta unas cuantas veces antes de rendirse, porque
# sin estas cuentas el paso 02c (islotes) no corre.
for a in GSE159984 GSE108539; do
  f=${a}_raw_counts_GRCh38_p13_NCBI.tsv.gz
  for intento in 1 2 3 4 5; do
    [ -s "$f" ] && gzip -t "$f" 2>/dev/null && break
    rm -f "$f"
    get "$N&acc=$a&format=file&file=${a}_raw_counts_GRCh38.p13_NCBI.tsv.gz" "$f"
    if [ -s "$f" ] && ! gzip -t "$f" 2>/dev/null; then
      rm -f "$f"
      echo "  [reintento $intento] $a: NCBI devolvio HTML (CAPTCHA)"
      sleep $((intento * 20))
    fi
  done
  if [ -s "$f" ] && gzip -t "$f" 2>/dev/null; then :; else
    echo "  [FALLO] $a: NCBI sigue devolviendo HTML. Reintente mas tarde, o para GSE159984 use"
    echo "          T2D_ARCHS4_DL=1 bash download_public.sh && python python/analyses/00c_archs4_counts.py"
  fi
done
[ -s GSE159984_raw_counts_GRCh38_p13_NCBI.tsv.gz ] && ln -sfn GSE159984_raw_counts_GRCh38_p13_NCBI.tsv.gz GSE159984_raw_counts.tsv.gz || true
get "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66446/suppl/GSE66446_Cushings_htseq_gene_counts_GRCh37.74.txt.gz" GSE66446_Cushings_htseq_gene_counts_GRCh37.74.txt.gz
# la serie tiene dos plataformas; 00_prepare_data.py busca el nombre sin sufijo
[ -s GSE159984-GPL16791_series_matrix.txt.gz ] && ln -sfn GSE159984-GPL16791_series_matrix.txt.gz GSE159984_series_matrix.txt.gz

echo "## 9. recount3: GSE108539 (SRP127581), sustituto de las cuentas NCBI de la plantilla de fase"
mkdir -p recount3
R=https://duffel.rail.bio/recount3/human/data_sources/sra
get $R/gene_sums/81/SRP127581/sra.gene_sums.SRP127581.G026.gz recount3/sra.gene_sums.SRP127581.G026.gz
get $R/metadata/81/SRP127581/sra.sra.SRP127581.MD.gz          recount3/sra.sra.SRP127581.MD.gz
get $R/metadata/81/SRP127581/sra.recount_qc.SRP127581.MD.gz recount3/sra.recount_qc.SRP127581.MD.gz
get $R/gene_sums/49/SRP055749/sra.gene_sums.SRP055749.G026.gz recount3/sra.gene_sums.SRP055749.G026.gz
get $R/metadata/49/SRP055749/sra.sra.SRP055749.MD.gz          recount3/sra.sra.SRP055749.MD.gz
get $R/metadata/49/SRP055749/sra.recount_qc.SRP055749.MD.gz recount3/sra.recount_qc.SRP055749.MD.gz

if [ -n "${T2D_ARCHS4_DL:-}" ]; then
  echo "## 10. ARCHS4 human_gene_v2 (62 GB), para GSE159984 (sin cuentas de autor ni en recount3)"
  mkdir -p archs4
  curl -s -S -L -C - -o archs4/human_gene_v2.latest.h5 https://s3.k8s.maayanlab.cloud/archs4/files/human_gene_v2.latest.h5 \
    && echo "  [ok]    archs4/human_gene_v2.latest.h5" || echo "  [FALLO] ARCHS4"
fi

# fuentes sin version (gene_info cambia a diario, OmniPath es una consulta en vivo): se registra fecha y sha256
{ echo "# $(date -u +%Y-%m-%dT%H:%MZ)"; sha256sum _soft/Homo_sapiens.gene_info.gz dorothea_ABC.tsv Human_GRCh38_p13_annot_tsv.gz 2>/dev/null; } >> CHECKSUMS.txt

echo "## 11. pyCHIRAL, para el paso 06 de concordancia de fase (no esta en PyPI)"
# Revision fijada: la que se uso para verificar el paso. PYCHIRAL_COMMIT permite cambiarla.
PYCHIRAL_COMMIT="${PYCHIRAL_COMMIT:-e1ea76e6f01c3106993fcf93b9195bb351904f58}"
if [ -d pychiral/.git ] && [ "$(git -C pychiral rev-parse HEAD 2>/dev/null)" = "$PYCHIRAL_COMMIT" ]; then
  echo "  [ya]    pychiral @ ${PYCHIRAL_COMMIT:0:7}"
else
  rm -rf pychiral
  if git clone -q https://github.com/AndreaSalati/pyCHIRAL.git pychiral 2>/dev/null \
     && git -C pychiral -c advice.detachedHead=false checkout -q "$PYCHIRAL_COMMIT" 2>/dev/null; then
    echo "  [ok]    pychiral @ ${PYCHIRAL_COMMIT:0:7}"
    echo "# pyCHIRAL $(git -C pychiral rev-parse HEAD)" >> CHECKSUMS.txt
  else
    rm -rf pychiral
    echo "  [FALLO] pyCHIRAL: el paso 06 dira que CHIRAL no se ejecuto"
  fi
fi

echo "### $(date)  fin"
echo "## Informe de 00_prepare_data.py:"
cd "$ROOT" && "${PY:-python3}" python/analyses/00_prepare_data.py
