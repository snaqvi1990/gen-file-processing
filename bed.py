#common functions for processing .bed files

#given a standard bed line, returns all info in a dictionary format. Hierarchical with the extra info at the end
def readline(bed_line,coords_as_ints = 'n'):
    add_info = {}
    add_fields = ['name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
    split_line = bed_line.strip().split('\t')
    chrom, chromStart, chromEnd = split_line[:3]
    if coords_as_ints == 'y':
        chromStart = int(chromStart)
        chromEnd = int(chromEnd)
    for i in range(len(split_line[3:])):
        add_info[add_fields[i]] = split_line[3:][i]
    line_info = {'chrom':chrom,'start':chromStart,'end':chromEnd,'add_info':add_info}
    return line_info

def writeline(line_info):
    add_fields = ['name','score','strand','thickStart','thickEnd','itemRgb','blockCount','blockSizes','blockStarts']
    string_towrite = []
    string_towrite.append(line_info['chrom'])
    string_towrite.append(line_info['start'])
    string_towrite.append(line_info['end'])
    for field in add_fields:
        if field in line_info['add_info']:
            string_towrite.append(line_info['add_info'][field])

    return '\t'.join(map(str,string_towrite)) + '\n'

def readfile(bed_file,list_or_dict = 'list',coords_as_ints = 'n'):
    if list_or_dict == 'list':
        bed = []
    elif list_or_dict == 'dict':
        bed = {}
    f = open(bed_file,'r')
    for line in f:
        line_info = readline(line,coords_as_ints)
        if list_or_dict == 'list':
            bed.append(line_info)
        elif list_or_dict == 'dict':
            bed[line_info] = 1
    f.close()
    return bed

# returns the length of overlap between two sequences given their coordinates as tuples (start,end), and also the coordinates of the overlap (0,0) if none
def overlap(seq1,seq2):
    start1 = int(seq1[0])
    end1 = int(seq1[1])
    start2 = int(seq2[0])
    end2 = int(seq2[1])
    if (start1 <= end2) & (start1 >= start2):
        overlap = abs(end2-start1)
        coords = (start1,end2)
    elif (start2 <= end1) & (start2 >= start1):
        overlap = abs(end1-start2)
        coords = (start2,end1)
    else:
        overlap = 0
        coords = (0,0)
    return overlap, coords

