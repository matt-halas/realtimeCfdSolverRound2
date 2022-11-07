NX = 40
NY = 40


def two_to_one(x_idx, y_idx):
    return x_idx + y_idx * NY

def one_to_two(idx):
    return idx % NX, idx // NY