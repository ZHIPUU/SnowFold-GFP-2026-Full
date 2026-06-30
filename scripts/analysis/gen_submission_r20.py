import json, csv

# Load R20 Top 6
with open(r'D:\workspace\round20\r20_top6.json') as f:
    top6 = json.load(f)

# Generate submission CSV
csv_path = r'D:\workspace\round20\submission_r20.csv'
with open(csv_path, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['Team_Name', 'Seq_ID', 'Sequence'])
    for i, c in enumerate(top6):
        w.writerow(['SnowFold', i+1, c['seq']])

print(f'Submission CSV: {csv_path}')
print(f'6 候选已写入')
