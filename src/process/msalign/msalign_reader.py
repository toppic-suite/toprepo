from torch.utils.data import Dataset as DataSet

class MsalignReader(DataSet):
    def __init__(self, msalign_file): 
        self.msalign_file = msalign_file

    def readmsalign_iter(self):
        f = open(self.msalign_file, "r", buffering=1024*1024*1024)
        current = None
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("BEGIN IONS"):
                current = {"meta": {}, "meta_lines": [], "peak_lines": []}
            elif line.startswith("END IONS"):
                if current:
                    yield current 
                    current = None
            elif current != None:
                if "=" in line:
                    key, val = line.split("=", 1)
                    current["meta"][key] = val
                    current["meta_lines"].append(line)
                else:
                    # Each line should be: mz intensity ion_type
                    current["peak_lines"].append(line)