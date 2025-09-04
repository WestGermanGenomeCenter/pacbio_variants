mkdir -p ~/gnomad_v4
cd ~/gnomad_v4
for CHR in {1..22} X Y; do
  curl -L -O "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.chr${CHR}.vcf.bgz"
  curl -L -O "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.chr${CHR}.vcf.bgz.tbi"
done