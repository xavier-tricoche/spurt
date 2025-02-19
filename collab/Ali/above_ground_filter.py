import numpy as np 
import nrrd 
import argparse
import tqdm


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Algorithms to filter out spurious above ground patterns in ground penetrating images')

    parser.add_argument('-i', '--input', type=str, required=True, help='Input NRRD image')
    parser.add_argument('-o', '--output', type=str, required=True, help='Output NRRD image')

    args = parser.parse_args()
    input, _ = nrrd.read(args.input)
    
    low = np.min(np.ravel(input))
    print(f'low={low}')
    
    w, h, d = input.shape

    elevation = np.argmax(input, axis=2)
    
    for n in tqdm.tqdm(range(0, w*h)):
        x = n // h
        y = n % h
        z = elevation[x,y]
        input[x,y,0:z] = low

    nrrd.write(args.output, input, compression_level=2)