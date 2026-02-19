class MsalignWriter():
    def __init__(self, msalign_file): 
        self.msalign_file = msalign_file
        self.f = open(msalign_file, "w", buffering=1024*1024*1024)
    
    def __del__(self):
        self.f.close()
    
    def close(self):
        self.f.close()

    def write(self, spectrum):
        self.f.write("BEGIN IONS\n")
        for line in spectrum["meta_lines"]:
            self.f.write(line + "\n")
        for line in spectrum["peak_lines"]:
            self.f.write(line + "\n")
        self.f.write("END IONS\n\n")

    def write_using_meta(self, spectrum):
        self.f.write("BEGIN IONS\n")
        for key in spectrum["meta"]:
            self.f.write(key + "=" + spectrum["meta"][key] + "\n")
        for line in spectrum["peak_lines"]:
            self.f.write(line + "\n")
        self.f.write("END IONS\n\n")
    def write_mz_intensity(self, spectrum):
        self.f.write("BEGIN IONS\n")
        for line in spectrum["meta_lines"]:
            self.f.write(line + "\n")
        for line in spectrum["peak_lines"]:
            fields = line.strip().split()
            mass = float(fields[0])
            charge = int(fields[2])
            mz = mass / charge + 1.007276466879
            # format mz with 5 fractional positions
            self.f.write(f"{mz:.5f}\t{fields[1]}\t{fields[2]}\t{fields[3]}\n")
        self.f.write("END IONS\n\n")