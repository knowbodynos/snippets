# s3

## Usage

```
FASTA_PATH = os.path.join(DATA_PATH, 'ref', 'Zea_mays.B73_RefGen_v4.fa')
GTF_PATH = os.path.join(DATA_PATH, 'ref', 'Zea_mays.B73_RefGen_v4.gtf')

S3_KEYS = {'maize_b73.agpv4/hisat2_spliced/Zea_mays.B73_RefGen_v4.dna_sm.toplevel.fa': FASTA_PATH,
           'maize_b73.agpv4/hisat2_spliced/Zea_mays.B73_RefGen_v4.59.gtf': GTF_PATH}

s3.pull(BUCKET_NAME, S3_KEYS)
```