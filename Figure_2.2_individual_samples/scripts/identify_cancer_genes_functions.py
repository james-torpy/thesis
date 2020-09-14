# for each DE gene, scrapes google scholar for articles:
def scholar_scrape(
    sample_id,
    genes,
    incl_terms,
    rm_terms,
    no_returned,
    out_dir
):
    import os
    from bs4 import BeautifulSoup
    import requests
    import re
    import time
    import sys
    import random
    
    # remove any existing DE gene article file:
    if os.path.isfile(out_dir + sample_id + '_DE_gene_articles.txt'):
        os.remove(out_dir + sample_id + '_DE_gene_articles.txt')
    # scrape google scholar for articles on genes containing key terms:
    print("Finding articles on DE genes from sample " + sample_id)
    
    for i in range(1, len(genes)):
        # determine waiting time at random, or wait 30 min every 20 entries:
        if i % 20 == 0:
            wait_time = 1800
        else:
            wait_time = random.randint(120,300)
        print('\nWaiting for ' + str(wait_time) + ' sec to avoid bot detection...\n')
        time.sleep(wait_time)
        print('\nScraping Google Scholar for articles on ' + genes[i] + \
            ' containing key search terms...')
        headers = {'User-Agent':'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_11_2) AppleWebKit/601.3.9 (KHTML, like Gecko) Version/9.0.2 Safari/601.3.9'}
        url = 'https://scholar.google.com.au/scholar?hl=en&as_sdt=0%2C5&q=' + \
            genes[i] + '&btnG='
        response=requests.get(url,headers=headers)
        soup=BeautifulSoup(response.content,'lxml')
    
        # store citations:
        cites = []
        for item in soup.select('[data-lid]'):
            try:
                #print('----------------------------------------')
                citn = item.select('h3')[0].get_text()
                #print(citn) 
                cites.append(citn)
            except Exception as e:
                #raise e
                print('')
        
        # remove anything in brackets from cites:
        cites = [re.sub('\[[^>]+\]', '', s) for s in cites]
        
        # add terms with capitals as first letter, and all capitals, to incl_terms:
        incl_terms.extend([t.upper() for t in incl_terms])
        incl_terms.extend([t.capitalize() for t in incl_terms])
        incl_terms = list(set(incl_terms))
        
        rm_terms.extend([t.upper() for t in rm_terms])
        rm_terms.extend([t.capitalize() for t in rm_terms])
        rm_terms = list(set(rm_terms))
        
        # search citations for key terms and add to list:
        key_cites = []
        for kterm in incl_terms:
            term_ind = [i for i, s in enumerate(cites) if kterm in s]
            if len(term_ind) > 0:
                key_cites.extend(list(map(cites.__getitem__, term_ind)))
        
        # remove citations containing bad keywords:
        key_cites = [x for x in key_cites if all(y not in x for y in rm_terms)]
        
        # keep only unique elements:
        key_cites = list(set(key_cites))

        print('\nCitations containing keywords are:\n')
        for c in key_cites:
            print(c.encode(sys.stdout.encoding, errors='replace'))
            print("\n")

        # if key term present, record n citations in text document:
        final_cites = '\n\n'.join(key_cites[0:no_returned])
        if len(final_cites) > 0:
            with open(
                out_dir + sample_id + '_DE_gene_articles.txt', mode='at', encoding='utf-8'
            ) as myfile:
                myfile.write('\n\n----------------------------------------')
                myfile.write('\n' + genes[i] + ':\n\n')
                myfile.write(final_cites)
    
            print('\n----------------------------------------\n')