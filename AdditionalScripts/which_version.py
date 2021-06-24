import os
import sys

def main():
    os.system(sys.argv[1]+ ' 2>&1')

if __name__ == '__main__':
    main()
