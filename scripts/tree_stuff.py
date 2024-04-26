import subprocess
import argparse


def tree(input, model):
    make_tree_file = f"iqtree2 -s {input} -m {model} -T AUTO"
    subprocess.run(make_tree_file, shell=True)


def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    parser.add_argument('-in', '--db_translated_file', help="db_translated_file", type=str)
    parser.add_argument('-m', '--model', help="MDF for modelfinder", type=str)
    args = parser.parse_args()

    db_translated_file = args.db_translated_file
    model = args.model

    tree(db_translated_file, model)


if __name__ == '__main__':
    main()
    