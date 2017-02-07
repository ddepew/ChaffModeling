# INTRO TO FILE
import math
import random
import numpy as np
from numpy import linalg as la



#################################################
# Read inputs from files
def main():

    inputs = open('input_positions.txt','r')

    simTag = inputs.readline().strip()
    print(simTag)

    numMol = int(inputs.readline().strip())
    print(numMol)

    # Initialize dictionary containing molecule info
    molTypes = {'Molecule_Name':[], 'Number_Molecules':[], 'Region_Occupation':[], 'Enforced_Boundaries':[]}

    for ii in range(0,numMol):
        molName = inputs.readline().strip()
        molTypes['Molecule_Name'].append(molName)

        molTypes['Number_Molecules'].append(int(inputs.readline().strip()))

        region = inputs.readline().strip().split(' ')
        molTypes['Region_Occupation'].append([float(bnd) for bnd in region])

        enforced = inputs.readline().strip().split(' ')
        molTypes['Enforced_Boundaries'].append([int(bnd) for bnd in enforced])

        # Would read in random or structured input here...
        inputs.readline()


    #print(molTypes)

    # Get global xmin from each molecule's region of occupation
    xmin = min([lst[0] for lst in molTypes['Region_Occupation']])
    xmax = max([lst[1] for lst in molTypes['Region_Occupation']])
    ymin = min([lst[2] for lst in molTypes['Region_Occupation']])
    ymax = max([lst[3] for lst in molTypes['Region_Occupation']])
    zmin = min([lst[4] for lst in molTypes['Region_Occupation']])
    zmax = max([lst[5] for lst in molTypes['Region_Occupation']])

    #print([xmin, xmax, ymin, ymax, zmin, zmax])


    # Open output files
    setup = open('setup.xyz','w')
    data = open('data.'+simTag,'w')

    # (uniAtoms, numAtoms) = atom(A)
    uniAtomsGlobal = []
    numAtomsGlobal = []

    for x in molTypes['Molecule_Name']:
        (uniAtoms,numAtoms) = atom(x)
        uniAtomsGlobal.append(uniAtoms)
        numAtomsGlobal.append(numAtoms)


   #calculate the total atoms (number of atoms per molecule*number of molecules) in the entire simulation
    simAtoms = sum([l*k for l,k in zip(numAtomsGlobal,molTypes['Number_Molecules'])])
    print(simAtoms)

    #calculate total number of unique atoms
    uniAtomsGlobal = sum(uniAtomsGlobal)

    #Header lines for setup
    print(simAtoms,file=setup)
    print(simTag,file=setup)

    #Header lines for data
    print('# ' + simTag, file=data, end='\n\n')
    print(str(simAtoms)+' atoms\n'+str(uniAtomsGlobal)+' atom types', file = data, end='\n\n')
    print(str(xmin)+' '+str(xmax)+' xlo xhi\n',file=data)
    print(str(ymin) + ' ' + str(ymax) + ' ylo yhi\n', file=data)
    print(str(zmin) + ' ' + str(zmax) + ' zlo zhi', file=data,end='\n\n')
    print('Masses',file=data,end='\n\n')

    for x in range(1,uniAtomsGlobal+1):
        print(x,file=data)

    print('\nAtoms',file=data,end='\n\n')
## Iterate through assigning positions for each molecule
    kk = 0
    for x in range(0,numMol):
        mol_xyz = [[],[],[]]
        mol_elem = []
       #Read in molecule base geometry
        try:
            molFile = open(molTypes['Molecule_Name'][x] +'.mol','r')

            for line in molFile:
                molCell = line.strip().split()
                print(molCell)
                mol_xyz[0].extend([float(molCell[0])])
                mol_xyz[1].extend([float(molCell[1])])
                mol_xyz[2].extend([float(molCell[2])])
                mol_elem.append(molCell[3])
        except FileNotFoundError:\

            mol_xyz = [[0.0,0.0,0.0]]
            mol_elem = [molTypes['Molecule_Name'][x]]
            print('error!')

        #Set boundaries for this atom
        xmin = molTypes['Region_Occupation'][x][0]
        xmax = molTypes['Region_Occupation'][x][1]
        ymin = molTypes['Region_Occupation'][x][2]
        ymax = molTypes['Region_Occupation'][x][3]
        zmin = molTypes['Region_Occupation'][x][4]
        zmax = molTypes['Region_Occupation'][x][5]

        restr = []

        for idim in range(0,2):
            Inf = float("inf")
            if molTypes["Enforced_Boundaries"][x][idim*2-1] == 1:
                restr.append([-Inf,Inf,-Inf,Inf,-Inf,Inf])
                restr[-1][idim*2] = molTypes["Enforced_Boundaries"][x][idim*2-1]

            if molTypes["Enforced_Boundaries"][x][idim*2] == 1:
                restr.append([-Inf,Inf,-Inf,Inf,-Inf,Inf])
                restr[-1][idim*2-1] = molTypes["Enforced_Boundaries"][x][idim*2]

        for imol in range(0,molTypes["Number_Molecules"][x]):
            pos = randomPosition(mol_xyz, xmin, xmax, ymin, ymax, zmin, zmax)
            bnd = [min([r[0] for r in pos]), max([r[0] for r in pos]),min([r[1] for r in pos]), max([r[1] for r in pos]),min([r[2] for r in pos]), max([r[2] for r in pos])]
            conf = False
            for sublist in restr:
                if overlap(bnd, sublist):
                    conf = True
            while conf:
                conf = False
                pos = randomPosition(mol_xyz, xmin, xmax, ymin, ymax, zmin, zmax)
                bnd = [min([r[0] for r in pos]), max([r[0] for r in pos]), min([r[1] for r in pos]),
                       max([r[1] for r in pos]), min([r[2] for r in pos]), max([r[2] for r in pos])]
                for sublist in restr:
                    if overlap(bnd, sublist):
                        conf = True
            restr.append(bnd)
            extraBnd = edgeBC(bnd, xmin, xmax, ymin, ymax, zmin, zmax, molTypes["Enforced_Boundaries"][x])

            if extraBnd:
                restr.append(extraBnd)

            #write to file and print to console
            for sub1,sub2 in zip(pos,mol_elem):
                kk +=  1
                print("%s %8.5f %8.5f %8.5f" %(sub2, sub1[0],sub1[1],sub1[2]), file=setup)
                print("%3d %2s 0.0 %8.5f %8.5f %8.5f" %(kk, sub2, sub1[0],sub1[1],sub1[2]),file=data)

            print('Atom positions for molecule' + str(imol) + ':')
            for  atm in pos:
                for el in atm:
                    print(str(el),end=' ')
                print('')
            print('')

def atom(A):

    uniAtoms = len([l for l in A if l.isupper()])
    B = [int(l) for l in A if l.isdigit()]
    totalAtoms = sum(B) + (uniAtoms - len(B))

    #totalAtoms = sum(int([l for l in A if l.isdigit()]))
    return (uniAtoms, totalAtoms)
# upper = []
# for letter in A:
#     if letter.isdigit():
#         upper.append(int(letter))
# sum(upper)

def matrix(R,v):
    rotated = []
    for sublist in R:
        rotated.append(sublist[0]*v[0] + sublist[1]*v[1] + sublist[2]*v[2])
    return rotated


def randomPosition(mol_xyz,xmin,xmax,ymin,ymax,zmin,zmax):

#Returns a random 3D position for a Hydrazine molecule with center of mass within the box defined by the input limits
    R = randRotation()
    rotated = []
    for atom in mol_xyz:
        rotated.append(matrix(R,atom))
    #generate a random orientation for the molecule
    #generate a random position for the origin (COM) of the molecule
    origin = [(xmax-xmin)*random.random()+xmin, (ymax-ymin)*random.random()+ymin, (zmax-zmin)*random.random()+zmin]

    mol = []
    for atom in rotated:
        mol.append([origin[0] + atom[0],origin[1] + atom[1],origin[2] + atom[2]])

    return(mol)

def randUnitVector():
    import math
    import random
    import numpy as np
    phi = 2*math.pi*random.random()
    theta = math.acos(1-2*random.random())
    u = np.array([(math.sin(theta)*math.cos(phi)),(math.sin(theta)*math.sin(phi)),(math.cos(theta))])

    return u
##Check that it is a unit vector
# u = randUnitVector()
# check = ((((u[0]) ** 2) + ((u[1]) ** 2) + ((u[2]) ** 2)) ** .5)
# print(u)
# print(check)


def cross(a, b):
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]

    return c


def randRotation():
    # startTime = timeit.default_timer()
    # print('Time after zero function:' + str(timeit.default_timer() - startTime))
    R = np.zeros((3,3),float)
    u1 = randUnitVector()
    R[:,0] = np.array([u1])
    u2 = randUnitVector()
    mu = u2 - (np.dot((np.dot(u2.transpose(), u1)),u1))

    while la.norm(mu) > 0.1:
        u1 = randUnitVector()
        R[:, 0] = np.array([u1])
        u2 = randUnitVector()
        mu = u2 - ((np.dot(u2.transpose(),u1)) * u1)

    R[:,1] = np.array([(mu / la.norm(mu))])
    R[:,2] = np.array([(np.cross(R[:,0],R[:,1]))])
    return(R)

def overlap(bnd1,bnd2):

    cond1 = bnd1[0] > bnd2[1]  #vol 1 left face is right of vol 2 right face
    cond2 = bnd2[0] > bnd1[1]  #vol 1 right face is left of vol 2 left face
    cond3 = bnd1[2] > bnd2[3]  #vol 1 front face is behind vol 2 back face
    cond4 = bnd2[2] > bnd1[3]  #vol 1 back face is in front of vol 2 front face
    cond5 = bnd1[4] > bnd2[5]  #vol 1 bottom face is above vol 2 top face
    cond6 = bnd2[4] > bnd1[5]  #vol 1 top face is below vol 2 bottom face

    b= not cond1 and not cond2 and not cond3 and not cond4 and not cond5 and not cond6
    return (b)

def edgeBC(bnd, xmin, xmax, ymin, ymax, zmin, zmax, bc):
    bnd_copy = bnd
    changed = False
    gblBnd = [xmin,xmax,ymin,ymax,zmin,zmax]

    for i in range(0,2):
        imin = 2*i
        imax = 2*i + 1
        if not bc[imin] and (bnd[imin] < gblBnd[imin]):
            bnd_copy[imin] = gblBnd[imin]
            bnd_copy[imax] = bnd[imax] - gblBnd[imax]
            changed = True
        if not bc[imax] and (bnd[imax] > gblBnd[imax]):
            bnd_copy[imin] = gblBnd[imin]
            bnd_copy[imax] = bnd[imax] - gblBnd[imax]
            changed = True

    if not changed:
        bnd_copy = False
    return(bnd_copy)


main()

