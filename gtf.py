#common functions for processing .gtf files

#given a standard gtf line, returns all info in a dictionary format. Hierarchical with the extra info at the end
def readline(gtf_line,coords_as_ints = 'n'):
    line_info = {}
    att_dict = {}
    att_order = []
    chrom, source, feature, start, end, score, strand, frame, attributes = gtf_line.strip().split('\t')
    attributes = attributes.strip(';')
    if coords_as_ints == 'y':
        start = int(start)
        end = int(end)
    att_pair_list = attributes.split('; ')
    if len(att_pair_list) == 1:
        att_dict = {'single':att_pair_list}
    else:
        for att_pair in att_pair_list:
            #print(att_pair.strip().split(' '))
            att, value = att_pair.strip().split(' \"')
            att_order.append(att)
            value = value.strip('\"')
            att_dict[att] = value
    line_info = {'chrom':chrom,'source':source,'feature':feature,'start':start,'end':end,'score':score,'strand':strand,'frame':frame,'attributes':att_dict,'attributes_order':att_order}
    return line_info


#given a dict of the info on a GTF line in the same way as read() outputs, returns a string of the line ready to write to file
def writeline(line_info):
    att_string = ''
    for att in line_info['attributes_order']:
        att_string += att + ' \"' + line_info['attributes'][att] + '\"; '
    str_to_print = '\t'.join([line_info['chrom'],line_info['source'],line_info['feature'],str(line_info['start']),str(line_info['end']),line_info['score'],line_info['strand'],line_info['frame'],att_string,'\n'])
    return str_to_print
