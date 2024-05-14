import numpy as np
import matplotlib.pyplot as plt
from Bio.Align import substitution_matrices

def nw(x, y, seq_type='DNA', match=1, mismatch=1, gap=2, sub_matrix=None):
    def populate_matrix(F, P):
        for i in range(1, nx + 1):
            for j in range(1, ny + 1):
                if seq_type == 'PROTEIN':
                    score = sub_matrix.get((x[i - 1], y[j - 1]), sub_matrix.get((y[j - 1], x[i - 1]), mismatch))
                else:  # DNA or RNA
                    if x[i - 1] == y[j - 1]:
                        score = match
                    else:
                        score = -mismatch
                
                t[0] = F[i - 1, j - 1] + score
                t[1] = F[i - 1, j] - gap
                t[2] = F[i, j - 1] - gap
                tmax = np.max(t)
                F[i, j] = tmax
                if t[0] == tmax:
                    P[i, j] += 2
                if t[1] == tmax:
                    P[i, j] += 3
                if t[2] == tmax:
                    P[i, j] += 4

    nx = len(x)
    ny = len(y)

    # Optimal score at each possible pair of characters.
    F = np.zeros((nx + 1, ny + 1))
    F[:, 0] = np.linspace(0, -nx * gap, nx + 1)
    F[0, :] = np.linspace(0, -ny * gap, ny + 1)

    # Pointers to trace through an optimal alignment.
    P = np.zeros((nx + 1, ny + 1))
    P[1:, 0] = 3
    P[0, 1:] = 4
    
    # Temporary scores.
    t = np.zeros(3)

    populate_matrix(F, P)

    # Trace through an optimal alignment.
    i = nx
    j = ny
    rx = []
    ry = []
    while i > 0 or j > 0:
        if P[i, j] == 2 or P[i, j] == 5 or P[i, j] == 6 or P[i, j] == 9:
            rx.append(x[i - 1])
            ry.append(y[j - 1])
            i -= 1
            j -= 1
        elif P[i, j] == 3 or P[i, j] == 5 or P[i, j] == 7 or P[i, j] == 9:
            rx.append(x[i - 1])
            ry.append('-')
            i -= 1
        elif P[i, j] == 4 or P[i, j] == 6 or P[i, j] == 7 or P[i, j] == 9:
            rx.append('-')
            ry.append(y[j - 1])
            j -= 1
        else:
            break

    # Reverse the strings.
    rx = ''.join(rx)[::-1]
    ry = ''.join(ry)[::-1]
    return F, '\n'.join([rx, ry]), rx, ry, P

def plot_matrix_with_arrows(F, x, y, P, aligned_x, aligned_y):
    plt.figure(figsize=(12, 8))

    # Plotting matrix
    plt.imshow(F, cmap='viridis', origin='upper')
    plt.colorbar(label='Score')
    plt.title('Needleman-Wunsch Matrix')
    plt.xlabel('Sequence Y')
    plt.ylabel('Sequence X')
    plt.xticks(np.arange(len(y) + 1), [''] + list(y))
    plt.yticks(np.arange(len(x) + 1), [''] + list(x))

    i, j = len(x), len(y)
    while i > 0 or j > 0:
        if P[i, j] in [2, 5, 6, 9]:  # Diagonal
            plt.arrow(j - 0.15, i - 0.15, -0.25, -0.25, head_width=0.05, head_length=0.05, fc='r', ec='r')
            plt.text(j, i, int(F[i, j]), color='white', ha='center', va='center')
            i -= 1
            j -= 1
        elif P[i, j] in [3, 5, 7, 9]:  # Vertical
            plt.arrow(j, i - 0.15, 0, -0.25, head_width=0.05, head_length=0.05, fc='r', ec='r')
            plt.text(j, i, int(F[i, j]), color='white', ha='center', va='center')
            i -= 1
        elif P[i, j] in [4, 6, 7, 9]:  # Horizontal
            plt.arrow(j - 0.15, i, -0.25, 0, head_width=0.05, head_length=0.05, fc='r', ec='r')
            plt.text(j, i, int(F[i, j]), color='white', ha='center', va='center')
            j -= 1

    plt.show()

def select_substitution_matrix(matrix_name):
    if matrix_name == "PAM1":
        return substitution_matrices.load('PAM1')  
    elif matrix_name == "PAM120":
        return substitution_matrices.load('PAM120')  
    elif matrix_name == "PAM250":
        return substitution_matrices.load('PAM250')  
    elif matrix_name == "B62":
        return substitution_matrices.load('BLOSUM62')  
    elif matrix_name == "B80":
        return substitution_matrices.load('BLOSUM80')  
    else:
        raise ValueError("Invalid substitution matrix name.")

# Input sequences and parameters
seq_type = input("Enter sequence type (DNA, RNA, or Protein): ").upper()
x = input("Enter sequence X: ").upper()
y = input("Enter sequence Y: ").upper()

if seq_type == 'PROTEIN':
    matrix_name = input("Enter substitution matrix name (PAM1, PAM120, PAM250, B62, B80): ").upper()
    sub_matrix = select_substitution_matrix(matrix_name)
else:
    sub_matrix = None

gap_score = abs(int(input("Enter gap score: ")))

matrix, alignment, aligned_x, aligned_y, pointers = nw(x, y, seq_type=seq_type, gap=gap_score, sub_matrix=sub_matrix)

print("Aligned Sequence X:", aligned_x)
print("Aligned Sequence Y:", aligned_y)

plot_matrix_with_arrows(matrix, x, y, pointers, aligned_x, aligned_y)
