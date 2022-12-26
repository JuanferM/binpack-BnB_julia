import os

pwd = os.getcwd()
filename = pwd.split('/')[-1]
instances = os.listdir()

dummy = 0
with open(filename+'.txt', 'w+') as output:
    output.write(str(len(instances)-1)+'\n')
    for ins in sorted(instances):
        if ins != 'script.py':
            output.write(' '+ins.split('.')[0]+" \n")
            with open(ins, 'r+') as f:
                lines = f.readlines()
                n = int(lines[0])
                C = int(lines[1])
                output.write(" ".join(['', str(C), str(n), str(dummy), '\n']))
                for i in range(2, len(lines)):
                    output.write(lines[i])
