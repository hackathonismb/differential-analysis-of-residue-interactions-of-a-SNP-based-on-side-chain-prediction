import csv
import json

if __name__ == "__main__":
    combined_ptms = {}
    with open('Combined_PTMs.csv', newline='') as csvfile:
        data = csv.reader(csvfile, delimiter=' ', quotechar='|')
        next(data)
        for row in data:
            uniprotID, _, PTMpos, PTMs, _ = row[0].rsplit(',')
            if uniprotID not in combined_ptms:
                combined_ptms[uniprotID] = {}
            if PTMs not in combined_ptms[uniprotID]:
                combined_ptms[uniprotID][PTMs] = []
            combined_ptms[uniprotID][PTMs].append(int(PTMpos))
    with open('Combined_PTMs.json', "w") as f:
        f.write(json.dumps(combined_ptms))
