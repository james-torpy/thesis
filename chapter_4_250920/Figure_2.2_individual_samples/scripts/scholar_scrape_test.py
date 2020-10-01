# -*- coding: utf-8 -*-
import os
from bs4 import BeautifulSoup
import requests
import re

headers = {'User-Agent':'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_11_2) AppleWebKit/601.3.9 (KHTML, like Gecko) Version/9.0.2 Safari/601.3.9'}
url = 'https://scholar.google.com.au/scholar?hl=en&as_sdt=0%2C5&q=mucl1&btnG='
response=requests.get(url,headers=headers)
soup=BeautifulSoup(response.content,'lxml')
print(soup.select('[data-lid]'))

## print and store first n citations:
#n_cites = []
#for item in soup.select('[data-lid]')[0:no_citations]:
#	try:
#		print('----------------------------------------')
#		citn = item.select('h3')[0].get_text()
#		print(citn)	
#		n_cites.append(citn)
#	except Exception as e:
#		#raise e
#		print('')
#print('----------------------------------------')
#
## remove anything in brackets from n_cites:
#n_cites = [re.sub('\[[^>]+\]', '', s) for s in n_cites]
#
## search citations for key terms and add to list:
#terms = []
#for kterm in key_terms:
#	term_ind = [i for i, s in enumerate(n_cites) if kterm in s]
#	if len(term_ind) > 0:
#		terms.extend(list(map(n_cites.__getitem__, term_ind)))
#
## keep only unique elements:
#terms = list(set(terms))
#print('\nSearch terms containing keywords are:\n')
#print(*terms, sep='\n\n')
#
## if key term present, record n citations in text document:
#with open(
#	out_dir + 'DE_gene_articles.txt', mode='wt', encoding='utf-8'
#) as myfile:
#    myfile.write('\n'.join(terms))
#






