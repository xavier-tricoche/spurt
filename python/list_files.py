import os
import argparse 

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Visualize path of particles')
    parser.add_argument('-i', '--input', type=str, required=True, help='Input filename')

    args = parser.parse_args()


    f = []
    for (dirpath, dirnames, filenames) in os.walk("//home//ed3brute//spurt//spurt-main//build"):
    #for (dirpath, dirnames, filenames) in os.walk("//mnt//c//Users//jacob//Desktop//lavd_results//reshaped - GAUSSIAN 1 3"):
        for file in filenames:
            if(args.input in file):
                f.append(file)
    f.sort()

    list_doc = open("//home//ed3brute//spurt//spurt-main//build//" + args.input + ".txt", "w")
    for line in f:
        list_doc.write(line + "\n")
    list_doc.close()
