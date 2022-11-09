NX = 40
NY = 40


def two_to_one(x_idx, y_idx):
    return x_idx + y_idx * NY

def one_to_two(idx):
    return idx % NX, idx // NY

def lerp(x1, x2, xp, y1, y2):
    return y1 + (y2 - y1) * (xp - x1) / (x2 - x1)