# forome_callers

Collection of variant callers based on BGM AB (allele balance) and bayesian callers.

This package includes a hraness to run BGM Rare Variant Callers. 
Currently implemented callers are:

**Allele Balance Callers**:

* De-Novo
* Homozygous Recessive
* Compound Heterozygous

**Bayesian Callers:**

* De-Novo 

The easiest way to run callers is to use variant_caller.py module.

```
usage: variant_caller.py [-h] -i VCF -f FAMILY [--results RESULTS]
                         [--dnlib DNLIB]

Run BGM variant callers

optional arguments:
  -h, --help            show this help message and exit
  -i VCF, --input VCF, --vcf VCF
                        Input VCF file, required. Use jointly called VCF for
                        better results
  -f FAMILY, --family FAMILY
                        Family (fam) file, required
  --results RESULTS     Results directory, only used for running Bayesian
                        callers
  --dnlib DNLIB         Path to De-Novo library. If specified, then Bayesian
                        De-Novo Caller is used, otherwise Allele Balance
                        Caller
```                        

Additional information on customized use of Bayesian De-Novo caller is in
package denovo2. In the same file is the instruction to build a custom    
De-Novo library. A prebuild library from prior BGM cases is available 
for download.                     

To include this variant caller in the pipeline see variant_caller.py 
and callers/harness.py modules
