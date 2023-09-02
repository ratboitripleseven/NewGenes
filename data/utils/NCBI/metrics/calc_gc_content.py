'''
calculate gc content

'''

def calc_gc_content(sequence)->float:
    g_count = 0
    a_count = 0
    c_count = 0
    t_count = 0
    
    for base in sequence:
        if base == 'G':
            g_count +=1
        elif base == 'A':
            a_count +=1
        elif base == 'C':
            c_count +=1
        elif base == 'T':
            t_count +=1
    
    
    return g_count, a_count, c_count, t_count, float(
        ((g_count+c_count)/(g_count+c_count+a_count+t_count))*100
    )
    
    

if __name__ == '__main__':
    g_c_content = calc_gc_content('ATGGAAAATATATTAGACCTGTGGAACCAAGCCC')
    print(g_c_content)