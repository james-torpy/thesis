# for each DE gene, scrapes google scholar for 
# articles:
def scholar_scrape(
    sample_id,
    genes,
    incl_terms,
    rm_terms,
    no_returned,
    out_dir
):
    # remove any existing DE gene article file:
    if os.path.isfile(out_dir + sample_name + '_DE_gene_articles.txt'):
        os.remove(out_dir + sample_name + '_DE_gene_articles.txt')
    #scrapes google scholar for articles on genes containing key terms:
    for de_gene in genes[0:3]:
        print("Finding articles on DE genes from sample " + sample_name)
        print('\nScraping Google Scholar for articles on ' + de_gene + \
            ' containing key search terms...')
        headers = {'User-Agent':'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_11_2) AppleWebKit/601.3.9 (KHTML, like Gecko) Version/9.0.2 Safari/601.3.9'}
        url = 'https://scholar.google.com.au/scholar?hl=en&as_sdt=0%2C5&q=' + \
            de_gene + '&btnG='
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
        #print(*key_cites.encode(sys.stdout.encoding, errors='replace'), sep='\n\n')
        
        # if key term present, record n citations in text document:
        if len(key_cites) > 0:
            with open(
                out_dir + sample_name + '_DE_gene_articles.txt', mode='at', encoding='utf-8'
            ) as myfile:
                myfile.write('----------------------------------------')
                myfile.write('\n' + de_gene + ':\n\n')
                myfile.write('\n\n'.join(key_cites[0:no_returned])
                myfile.write('\n\n')
            )
        
            print('\nWaiting for 45 sec to avoid bot detection...')
            time.sleep(45)
    
            print('\n----------------------------------------\n')