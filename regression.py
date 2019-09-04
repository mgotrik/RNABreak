import os


os.system("python partial.py > new_output.txt")

if os.path.exists("gold_output.txt"):
    print("Comparing...")
    os.system("diff gold_output.txt new_output.txt")
else:
    print("Reference does not exist; generating it.")
    os.system("mv new_output.txt gold_output.txt")
