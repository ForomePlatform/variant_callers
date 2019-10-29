from typing import Dict, Set, Tuple, List
import vcf as pyvcf


class AbstractCaller():
    def __init__(self):
        self.family = None
        self.samples = None
        self.unrelated_samples = set()
        return

    def init(self, family: Dict, samples: Set):
        self.family = family
        self.samples = samples
        self.unrelated_samples = self.samples - self.family.keys()
        self.affected_samples = [id for id in self.family if self.family[id]['affected']]
        self.unaffected_samples = [id for id in self.family if not self.family[id]['affected']]
        return

    def make_call(self, record: pyvcf.Reader) -> Tuple:
        return ()

    def get_tag(self):
        return __class__

    def get_type(self):
        return "String"

    def get_n(self):
        type = self.get_type()
        if type == "Flag":
            return 0
        return 1

    def get_description(self):
        return ""

    def get_header(self):
        pattern = '##INFO=<ID={tag},Number={n},Type={type},Description="{desc}">'
        return pattern.format(tag=self.get_tag(), type=self.get_type(),
                              n=self.get_n(), desc=self.get_description())

    @staticmethod
    def gt_to_int(gl: List):
        return [g if g != None else 0 for g in gl]