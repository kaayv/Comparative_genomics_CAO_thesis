import json
from tabulate import tabulate

def extract_info(filename, search_terms):
    init_annot = []
    with open(filename) as f:
        for jsonobj in f:
            #jsonobj1 = jsonobj.replace("\'", "\"")
            eggnog_dict = json.loads(jsonobj)
            init_annot.append(eggnog_dict)

#hold
    result = []
    headers = []
    for item in eggnog_dict:
        row = []
        for key, value in item.items():
            for search_term in search_terms:
                if search_term in str(value):
                    if not headers:
                        headers = ['Key'] + search_terms
                    row.append(value)
        if row:
            result.append([item['Key']] + row)

    print(tabulate(result, headers=headers))
#hold
def main():
    filename = 'example.json'
    search_terms = ['toppings', 'size']
    extract_info(filename, search_terms)

if __name__ == '__main__':
    main()
