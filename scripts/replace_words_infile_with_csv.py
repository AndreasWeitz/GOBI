import csv
import argparse
import re


def replace_words(csv_file, input_file, output_file):
    with open(csv_file, 'r', newline='', encoding='utf-8') as csv_input:
        reader = csv.reader(csv_input)
        word_dict = {row[0]: row[1] for row in reader}

    def replace_word(match, replacement):
        if match:
            return replacement + "  " + match
        else:
            return replacement

    with open(input_file, 'r', encoding='utf-8') as file:
        data = file.read()
        for word, replacement in word_dict.items():
            if re.match(r'^[A-Z]{7}', word):
                pattern = re.compile(r'\b' + word[:7] + r'(\w*)\b', re.IGNORECASE)
            else:
                pattern = re.compile(r'\b' + re.escape(word) + r'\b', re.IGNORECASE)
            data = pattern.sub(lambda match, replacement=replacement: replace_word(match.group(0), replacement), data)

    with open(output_file, 'w', encoding='utf-8') as modified_file:
        modified_file.write(data)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
    parser.add_argument('-csv', '--csv', help="Path to input csv file", type=str)
    parser.add_argument('-in', '--input_file', help="Path to file to change names", type=str)
    parser.add_argument('-out', '--output_file', help="output file with changed names", type=str)
    args = parser.parse_args()

    csv_file = args.csv
    input_file = args.input_file
    output_file = args.output_file

    replace_words(csv_file, input_file, output_file)
