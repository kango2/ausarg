while IFS=$'\t' read -r speciesid speciesname txid; do
    filename="${outdir}/${speciesid}.${speciesname// /_}.txid${txid}"
    echo datasets download genome accession ${speciesid} --include genome --filename ${filename}
done < /g/data/xl04/genomeprojects/referencedata/squamata.20241209.test.tsv
