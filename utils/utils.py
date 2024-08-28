def generate_path(disp_level: list, n: int=200):
    u = []
    for i, disp in enumerate(disp_level[1:]):
        for j in range(n):
            ui = disp_level[i] + (disp_level[i + 1] - disp_level[i]) * j / n
            u.append(ui)
    else:
        u.append(disp_level[-1])
    u = [round(i, 6) for i in u]
    return u

