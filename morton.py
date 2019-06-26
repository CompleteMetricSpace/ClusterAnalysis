#See wikipeida for more information

#Given a vector in [-1,1]^d, compute its morton code
def morton_code(vector):
    d = vector.__len__()
    vector = tuple((vector[i]+1)/2 for i in range(d))
    code = ''
    current_vector = list(vector)
    for i in range(32):
        for j in range(d):
            current_vector[j] = 2*current_vector[j]
            a = current_vector[j] // 1
            if a > 0:
                code += '1'
            else:
                code += '0'
            current_vector[j] = current_vector[j] - a
    return int(code, 2)
