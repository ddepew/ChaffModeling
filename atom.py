def atom(A):

    uniAtoms = len([l for l in A if l.isupper()])
    totalAtoms = sum([int(l) for l in A if l.isdigit()])
    #totalAtoms = sum(int([l for l in A if l.isdigit()]))
    return (uniAtoms, totalAtoms)

A = 'N2H4'
(uniAtoms, totalAtoms) = atom(A)
print(uniAtoms,totalAtoms)

# upper = []
# for letter in A:
#     if letter.isdigit():
#         upper.append(int(letter))
# sum(upper)

periodic_table = {}
