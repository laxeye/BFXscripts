#!/usr/bin/env python3
import requests
import xml.etree.ElementTree as ET
import ftplib
import sys

def connect_ftp():
    FTP_HOST = 'ftp.ncbi.nlm.nih.gov'
    FTP_USER = 'anonymous'
    FTP_PASS = 'anonymous@mail.com'

    ftp = ftplib.FTP(FTP_HOST)
    ftp.login(FTP_USER, FTP_PASS)

    return ftp

acc_file = sys.argv[1]
url_base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'

with open(acc_file) as acc_fh:
    acc_list = [x.rstrip() for x in acc_fh.readlines()]
    i = 0
    ftp = connect_ftp()
    for acc in acc_list:
        i += 1
        if i % 25 == 0:
            ftp.quit()
            ftp = connect_ftp()
        url = f'{url_base}esearch.fcgi?db=assembly&term={acc}&usehistory=y'
        with requests.get(url) as response:
            tree = ET.fromstring(response.text)
            key = tree[3].text
            web = tree[4].text
            #id = tree[5][0].text
            
            url2 = f'{url_base}esummary.fcgi?db=assembly&query_key={key}&WebEnv={web}'
            summary = ET.fromstring(requests.get(url2).text)
            ftp_folder = summary[0][1][46].text.split('/', maxsplit=3)[3]
            name = summary[0][1][6].text
            tax = summary[0][1][10].text
            f_name = f'{acc}_{name}_genomic.fna.gz'
            #ftp_path = f'{ftp_folder}/{f_name}'
            ftp.cwd('/')
            print(ftp_folder, summary[0][1][10].text)
            ftp.cwd(ftp_folder)
            with open(f_name, "wb") as file:
                try:
                    ftp.retrbinary(f"RETR {f_name}", file.write)
                except Exception as e:
                    print(f'File not found for accession\t{acc}')

ftp.quit()