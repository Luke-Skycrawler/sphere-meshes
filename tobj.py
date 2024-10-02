import numpy as np

def export_tobj(file, V, T):
    to_strv = lambda x: f"v {x[0]} {x[1]} {x[2]}\n"
    to_strt = lambda t: f"t {t[0]} {t[1]} {t[2]} {t[3]}\n"

    linesv = [to_strv(v) for v in V] 
    linest = [to_strt(t) for t in T]
    with open(file, 'w') as f:
        f.writelines(linesv + linest)

        