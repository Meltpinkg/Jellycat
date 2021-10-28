class Record(object):
    def __init__(self, record, idx):
        self.source = idx
        self.name = record.id
        self.type = 'BND' if record.info['SVTYPE'] == 'TRA' else record.info['SVTYPE']
        self.start = parse_to_int(record.pos)
        if self.type not in ['INS','DEL','INV','DUP','BND']:
            if 'INS' in record.alts[0]:
                self.type = 'INS'
            elif 'DEL' in record.alts[0]:
                self.type = 'DEL'
        if self.type == 'BND':
            if 'END' in record.info:
                self.end = parse_to_int(record.info['END'])
            else:
                try:
                    self.end = parse_to_int(record.stop)
                except:
                    self.end = 0
        if self.type != 'BND':
            if 'SVLEN' in record.info:
                self.end = abs(parse_to_int(record.info['SVLEN']))
            else:
                try:
                    self.end = abs(parse_to_int(record.stop) - self.start)
                except:
                    self.end = 0
        if self.type == 'BND':
            if 'CHR2' in record.info:
                chrom2 = record.info['CHR2']
            if 'END' in record.info:
                end = parse_to_int(record.info['END'])
            try:
                alt = str(record.alts[0])
                if alt[0] == ']':
                    tra_type = ']]N'
                    chrom2 = alt.split(':')[0][1:]
                    end = int(alt.split(':')[1][:-2])
                elif alt[0] == '[':
                    tra_type = '[[N'
                    chrom2 = alt.split(':')[0][1:]
                    end = int(alt.split(':')[1][:-2])
                else:
                    # print(type(alt))
                    if alt[1] == ']':
                        tra_type = 'N]]'
                        chrom2 = alt.split(':')[0][2:]
                        end = int(alt.split(':')[1][:-1])
                    else:
                        tra_type = 'N[['
                        chrom2 = alt.split(':')[0][2:]
                        end = int(alt.split(':')[1][:-1])
            except:
                print(record)
            self.chrom2 = chrom2
            self.end = end
            self.strand = tra_type
        '''
        if self.type == 'BND':
            tra_alt = str(record.alts[0])
            if tra_alt[0] == 'N':
                if tra_alt[1] == '[':
                    tra_type = 'A'
                else:
                    tra_type = 'B'
            elif tra_alt[0] == '[':
                tra_type = 'C'
            else:
                tra_type = 'D'
            if tra_alt[0] == 'N':
                tra_alt = tra_alt[2:-1]
            else:
                tra_alt = tra_alt[1:-2]
            self.chrom2 = tra_alt.split(':')[0]
            self.end = int(tra_alt.split(':')[1])
            self.strand = tra_type
        '''
        self.chrom1 = record.chrom
        self.strand = '.'
        if self.type != 'BND':
            self.chrom2 = record.chrom
            if 'STRAND' in record.info:
                self.strand = record.info['STRAND']
            elif 'STRANDS' in record.info:
                self.strand = record.info['STRANDS']
        if isinstance(self.strand, list) or isinstance(self.strand, tuple):
            self.strand = str(self.strand[0])
        self.ref = record.ref
        self.alt = record.alts  # tuple
        self.qual = record.qual  # NoneType
        if 'GT' in record.format:
            if record.samples[0]['GT'] == None or record.samples[0]['GT'][0] == None:
                self.gt = './.'
            else:
                self.gt = str(record.samples[0]['GT'][0]) + '/' + str(record.samples[0]['GT'][1])
        else:
            self.gt = './.'
        #if self.type == 'DEL':
        #    self.end = self.start + self.end

    def to_string(self):
        return self.name + ', ' + self.type + ', start: ' + str(self.start) + ', end: ' + str(self.end) + ', strand: ' + self.strand + ', gt: ' + self.gt + ', source: ' + str(self.source)

    def __eq__(self, other):
        if isinstance(other, Record):
            #if self.chrom == other.chrom and self.start == other.start and self.end == other.end and self.chrom2 == other.chrom2:
            if self.start == other.start and self.end == other.end:
                return True
            else:
                return False
        else:
            raise Exception("The type of object is not Record")
    def __lt__(self, other):
        if isinstance(other, Record):
            if self.start < other.start or (self.start == other.start and self.end < other.end):
                return True
            else:
                return False
        else:
            raise Exception("The type of object is not Record")
    def __gt__(self, other):
        if isinstance(other, Record):
            if self.start > other.start or (self.start == other.start and self.end > other.end):
                return True
            else:
                return False
        else:
            raise Exception("The type of object is not Record")

def parse_to_int(sth):
    if sth == None:
        return 0
    elif isinstance(sth, str):
        return int(sth)
    elif isinstance(sth, list):
        return parse_to_int(sth[0])
    elif isinstance(sth, tuple):
        return parse_to_int(sth[0])
    elif isinstance(sth, int):
        return sth
    else:
        return sth

def parse_svtype(sv_type):
    if 'DEL' in sv_type:
        return 'DEL'
    elif 'INS' in sv_type:
        return 'INS'
    elif 'INV' in sv_type:
        return 'INV'
    elif 'DUP' in sv_type:
        return 'DUP'
    elif 'BND' in sv_type or 'TRA' in sv_type:
        return 'BND'
    else:
        return sv_type