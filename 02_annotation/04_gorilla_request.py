import os
import sys
import subprocess
from selenium import webdriver
import pyautogui
import time

output_path = "/path/to/project/GORILLA-out"
os.chdir(output_path)
samples = os.listdir("/path/to/project/snpeff-out")

for s in samples:

    print(s)

    driver = webdriver.Firefox()
    driver.get('http://cbl-gorilla.cs.technion.ac.il/')

    upload = driver.find_element_by_name("target_file_name").send_keys("/path/to/project/snpeff-out/"+s+"/gene_list_extract_10")

    button = driver.find_element_by_name("run_gogo_button")
    button.click()

    time.sleep(5)

    pyautogui.hotkey('ctrl', 's')
    time.sleep(1)
    pyautogui.typewrite(s+"_10")
    pyautogui.hotkey('enter')

    time.sleep(3)
    os.system('mv /home/alex/Downloads/'+s+'_10.html '+output_path)
    os.system('mv /home/alex/Downloads/'+s+'_10_files/ '+output_path)
    #with open(output_path+s+".html", "w") as f:
        #f.write(driver.page_source)



















"""
import sys
import os
import requests
#import webbrowser

url = 'http://cbl-gorilla.cs.technion.ac.il/'
file = {'target_file_name': open('/path/to/project/snpeff-out/DB145_DB146/gene_list_extract', 'rb')}
submit = {'run_go_button': 'TRUE'}

x = requests.get(url, files = file, data = submit)

print(x.status_code)
print(x.text)
"""
