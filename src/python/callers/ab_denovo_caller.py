from typing import Dict, Set, Tuple
import vcf as pyvcf

from .ab_caller import ABCaller


class ABDenovoCaller(ABCaller):
    def __init__(self, recall_genotypes = True):
        super().__init__(recall_genotypes)

    def make_call(self, record: pyvcf.Reader) -> Tuple:
        genotypes = self.get_genotypes(record)
        if (not self.validate(genotypes)):
            return ()
        a = self.affected(genotypes)
        u = self.unaffected(genotypes)
        if (all([g > 0 for g in a]) and sum(u) == 0):
            return (self.get_tag(), None)
        return ()

    def get_tag(self):
        return "BGM_DE_NOVO"

    def get_type(self):
        return "Flag"

    def get_description(self):
        return "De novo by BGM allele balance caller"
