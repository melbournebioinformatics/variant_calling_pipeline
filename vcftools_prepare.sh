bgzip -c $1 > $1".gz"
tabix -p vcf $1".gz"
