## USAGE:

```sh
gff column types: seqid, source, type, pos, endpos, score, strand, phase, attr

gff3anno-linux-x86_64
-gff [path/to/file.gff3] #input gff3 file
-type {bed,vcf} #input file type (default: get from extension)
-skip N #(optional) for bed-file only, skip lines (default 0)
-header N #for bed-file only, header line num (default -1, no header)
-threads N #max threads (default 4)
-in path/to/input.{bed,vcf} #input bed or vcf file (or '-' for stdin, '-type' required)
-out path/to/output.{bed,vcf} #output file path or '-' for stdout
-seqid N #sequence id column number (default 1)
-pos N #position column number (default 2)
-endpos N #(optional) end position column number (info:<name> for vcf)
-where <par1>...<parN> #(optional) select from gff parameter (format <coltype>[:<attrname>]:<value>
-add <par1>...<parN> # fields to add to output file (format: <coltype>[:<attrname>])
```

## USAGE EXAMPLE:

```sh
gff3anno-linux-x86_64 -h
gff3anno-linux-x86_64 -gff annotation.gff3.gz -in test1.bed -out - -seqid 1 -pos 2 -where type:gene attr:gene_name:ADA -add attr:gene_name attr:gene_id type
```

## USAGE:

### add gene_name and gene_id to bed file from gencode gff3 file

```sh
$ cat test1.bed
#CHROM  START   END
chr1    10468   297829

$ gff3anno -gff gencode.v47.primary_assembly.basic.annotation.gff3.gz -endpos 3 -header 1 -threads $(nproc) -in test1.bed -out - -where type:gene -add attr:gene_name attr:gene_id
#CHROM  START   END     attr:gene_name  attr:gene_id
chr1    10468   297829  DDX11L16,DDX11L1,WASH7P,WASH7P,MIR6859-1,MIR1302-2HG,MIR1302-2,FAM138A,ENSG00000308361,ENSG00000290826,OR4G4P,OR4G11P,OR4F5,ENSG00000241860,ENSG00000239945,ENSG00000308314,CICP27,ENSG00000308579,ENSG00000268903,ENSG00000269981,ENSG00000239906,ENSG00000310528,RNU6-1100P,ENSG00000241599,DDX11L2,WASH9P,DDX11L17,WASH9P,MIR6859-2,ENSG00000308625,ENSG00000308544,ENSG00000292994,ENSG00000228463,ENSG00000306775,ENSG00000286448    ENSG00000290825.2,ENSG00000223972.6,ENSG00000310526.1,ENSG00000227232.6,ENSG00000278267.1,ENSG00000243485.6,ENSG00000284332.1,ENSG00000237613.3,ENSG00000308361.1,ENSG00000290826.2,ENSG00000268020.3,ENSG00000240361.3,ENSG00000186092.7,ENSG00000241860.8,ENSG00000239945.1,ENSG00000308314.1,ENSG00000233750.3,ENSG00000308579.1,ENSG00000268903.1,ENSG00000269981.1,ENSG00000239906.1,ENSG00000310528.1,ENSG00000222623.1,ENSG00000241599.1,ENSG00000308415.1,ENSG00000310527.1,ENSG00000279928.2,ENSG00000279457.4,ENSG00000273874.1,ENSG00000308625.1,ENSG00000308544.1,ENSG00000292994.2,ENSG00000228463.11,ENSG00000306775.1,ENSG00000286448.2
```

