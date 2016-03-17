###########################
### Functions: Methods: ###

# setting up with your email
class setup():
    # setup the email
    def email(self):
        email = input("Please enter your email: ")
        print ("You enter: ", email)
        return(email)
    # returns list of keywords you are interested in
    def keyword(self, numbers):
        keyword_lst = []
        while len(keyword_lst) != numbers:
            keyword = raw_input("Input keyword: ")
            keyword_lst.append(keyword)
        return(keyword_lst)

# retrieve the data
class retrieve():
    # returns list of keywords you are interested in
    def keyword(self, numbers):
        keyword_lst = []
        while len(keyword_lst) != numbers:
            keyword = raw_input("Input keyword: ")
            keyword_lst.append(keyword)
        return(keyword_lst)

# splitting a list in half iterable
def split_every(n, iterable):
    i = iter(iterable)
    piece = list(islice(i, n))
    while piece:
        yield piece
        piece = list(islice(i, n))

# list of records from the splitted pmid ids: this is the right one -- use parallel computing for this
def list_records(split_pmid):
    record_list = []
    for splits in range(len(split_pmid)):
        list_pmid = split_pmid[splits]
        temp = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs",from_uid=list_pmid))
        record_list.append(temp)

    return(record_list)

# Getting data from Medline: journals, authors, location
def record_fetch(pmid):
    # timer
    start_time_1 = timeit.default_timer()

    # extract data
    handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline",retmode="text")

    # timer
    print "The handle program time is ", timeit.default_timer() - start_time_1

    # use Medline parser to turn data into a parser object
    records = Medline.parse(handle)
    records = list(records)

    # timer
    start_time_2 = timeit.default_timer()

    # init data columns
    record_authors = []
    record_first = []
    record_journals = []
    record_dp = []
    record_place = []
    record_mesh = []
    record_pmid = []
    record_aff = []

    # iterate over records
    for record in records:
        record_authors.append(record.get("AU"))
        record_first.append(record.get("AU")[0])
        record_journals.append(record.get("JT"))
        record_dp.append(record.get("DP"))
        record_place.append(record.get("PL"))
        record_mesh.append(record.get("MH"))
        record_pmid.append(record.get("PMID"))
        record_aff.append(record.get("AD"))

    # timer
    print "The author program time is ", timeit.default_timer() - start_time_2

    # data frame
    data = pd.DataFrame({'authors': record_authors, 'first author': record_first, 'journals': record_journals, 'date publication': record_dp,
    'place': record_place, 'mesh term': record_mesh, 'pmid': record_pmid,'affliation': record_aff}, columns=['authors', 'first author', 'journals',
    'date publication', 'place', 'mesh term', 'pmid', 'affliation'])

    return(data)

# Getting data from Medline: first_author
def record_fetch_investigator(pmid_good):
    handle = Entrez.efetch(db="pubmed", id=pmid_good, rettype="medline",retmode="text")
    records = Medline.parse(handle)
    records = list(records)

    authors = []
    #first_author = []

    for record in records:
        authors.append(record.get("AU"))

    first_author = []

    for lst in range(len(authors)):
        if authors[lst] != None:
            index = authors[lst]
            first = index[0]
            first_author.append(first)

    return(first_author)

# Checks for connections in the record (just a list is passed)
def connections(record):
    # Cleaning out journals with none and keeping track of index
    record_clean= []
    counter_good_x = []
    counter_bad_x = []
    for x in range(len(record)):
        if record[x] != None:
            record_clean.append(record[x])
            counter_good_x.append(x)
        else:
            counter_bad_x.append(x)
    # Checking if there are any duplicates
    there_are_duplicates = len(record_clean)!=len(set(record_clean))
    # Check which authors are duplicated (show up more than once)
    record_duplicate_len = len([k for k,v in Counter(record_clean).items() if v>1])
    # To know the indices of duplicate authors
    D = defaultdict(list)
    for x,item in enumerate(record_clean):
        D[item].append(x)
    D = {k:v for k,v in D.items() if len(v)>1}

    return(there_are_duplicates, record_duplicate_len, D, counter_good_x, counter_bad_x, record_clean)

# Checks for author connections from pmid_good (uses list of lists)
def author_connections(record_authors):
    # Cleaning out authors with none and keeping track of index
    author_good = []
    counter_good_x = []
    counter_bad_x = []
    for x in range(len(record_authors)):
        if record_authors[x] != None:
            author_good.append(record_authors[x])
            counter_good_x.append(x)
        else:
            counter_bad_x.append(x)

    # Flattening the author list
    author_flat = list(itertools.chain.from_iterable(author_good))
    # Checking if there are any duplicates
    there_are_duplicates = len(author_flat)!=len(set(author_flat))
    # Check which authors are duplicated (show up more than once)
    author_duplicate_len = len([k for k,v in Counter(author_flat).items() if v>1])
    # To know the indices of duplicate authors
    D = defaultdict(list)
    for x,item in enumerate(author_flat):
        D[item].append(x)
    D = {k:v for k,v in D.items() if len(v)>1}

    return(there_are_duplicates, author_duplicate_len, D, counter_good_x, counter_bad_x, author_good)

# Make an edgelist (a list) I don't think you need to this for a just a list you can just use dict(zip ...)(Max)
def create_edgelist(investigators, article_id):

    combined = zip(investigators, article_id)

    investigator_to_articles = {}
    for combo in combined:
        investigator = combo[0]
        article = combo[1]

        if investigator not in investigator_to_articles:
            investigator_to_articles[investigator] = [article] #article_id goes in there
        else:
            investigator_to_articles[investigator].append(article)

    return(investigator_to_articles)

# Make a weighted graph (for a list) (max) not working right now
def create_weighted_graph(articles):
    G = nx.Graph()
    for article in range(len(articles)):
        new_edges = list(itertools.combinations(articles[article], 2))
        for new_e in new_edges:
            if not G.has_edge(new_e[0], new_e[1]):
                G.add_edge(new_e[0], new_e[1], weight=1)
            else:
                G[new_e[0]][new_e[1]]["weight"] += 1
        G.edges(data=True)

    return(G)

# Initial pmc_ids_temp list for each list that you want: doesn't work because records should be len of split not len of pmid
def pmc_ids(split_pmid, records, number):
    pmc_ids = []
    for record in range(len(records)):
        pmid_temp = split_pmid[record]
        records_temp = records[record]
        pmc_ids_temp = good_vs_bad(pmid_temp, records_temp, number)[0]
        pmc_ids.append(pmc_ids_temp)

    return(pmc_ids)

# get the pmc_ids and counter iteratively
def pmc_and_counter(split_pmid, records, number):
    good_vs_bad_list = []
    counter_good = []
    for record in range(len(records)):
        pmc = good_vs_bad(split_pmid[record], records[record], number)[0] #function call
        if len(pmc) != 0:
            good_vs_bad_list.append(pmc)
        counter = good_vs_bad(split_pmid[record], records[record], number)[1] #function call
        counter_good.append(counter)

    return(good_vs_bad_list, counter_good)

# get the pmid_good
def pmid_good_iter(counter_good, split_pmid):
    pmid_good = []
    for split in range(len(split_pmid)):
        temp = pmid_gooder(counter_good[split], split_pmid[split])
        pmid_good.append(temp)

    pmid_flat = list(itertools.chain.from_iterable(pmid_good))

    return(pmid_flat)

# manual way to save the articles collected from the scrap
def save_articles(articles):
    with open('keyword.txt','w') as file:
        for item in articles:
            print>>file, item + ","

# implementing complexity into the algorithm
def complexity(object):
    start_time = timeit.default_timer()
    object.iterator() # fix this
    print "The iterator program takes this much time: ", timeit.default_timer() - start_time
    start_time = timeit.default_timer()
    object.search() # fix this
    print "The search program takes this much time: ", timeit.default_timer() - start_time

# run analysis and create a list of suitable article ids - this code takes way too long - why? much longer than when i execute indiv
def analysis(keywords):
    article_lst = []
    ############################
    # [ initial handler setup ]
    for word in range(len(keywords)):
        print(keywords[word])
        # to find keyword population with count
        start_time = timeit.default_timer()
        handle = Entrez.esearch(db="pubmed",term=keywords[word],retmax=1)
        #handle = Entrez.esearch(db="pubmed",term=keywords[0], retmax=1)
        #handle = Entrez.esearch(db="pubmed",term=keywords[1],retmax=1)
        record_term = Entrez.read(handle)
        count = record_term['Count']

        if int(count) <= 10000:
            sample = 10000
        else:
            sample = float(count) * 0.30
            sample= int(sample)

        # taking 10% of population and using this as sample
        handle = Entrez.esearch(db="pubmed",term=keywords[word],retmax=sample)
        print "The initial handle program time: ", timeit.default_timer() - start_time

        ############################
        # [ handle test data ]
        #handle = Entrez.esearch(db="pubmed",term=keywords[0],retmax=sample)
        #handle = Entrez.esearch(db="pubmed",term=keywords[1],retmax=sample)

        print(sample)

        # storing the sample as the idlist for the articles

        ############################
        # [ splitting step ]
        start_time = timeit.default_timer()
        record_term = Entrez.read(handle)
        term_data = record_term["IdList"]
        pmid = term_data
        pmid = sorted(pmid)

        # splitting the record of articles for easier computation
        split_pmid = list(split_every(100, pmid))
        records = list_records(split_pmid)
        print "The split program time ", timeit.default_timer() - start_time

        ############################
        # [ data cleaning step ]
        # collecting citation data (pmc_ids) and cleaning
        start_time = timeit.default_timer()
        pmc_ids = pmc_and_counter(split_pmid, records, 5)[0]
        counter_good = pmc_and_counter(split_pmid, records, 5)[1]
        pmid_good = pmid_good_iter(counter_good, split_pmid)
        print "The data cleaning step program time: ", timeit.default_timer() - start_time

        ############################
        # [ fetch step ]
        # fetching author data and investigator data
        start_time = timeit.default_timer()
        #record_authors = record_fetch(pmid_good)[0]
        #authors = author_connections(record_authors)[5]

        # investigator (has be executed after authors)
        record_investigators = record_fetch_investigator(pmid_good)
        investigator_countg = connections(record_investigators)[3]
        investigators = connections(record_investigators)[5]

        # cleaning to get the true pmid_good after removing the none values in author, journal, and location
        article_id = pubmed_gooder(pmid_good, investigator_countg)
        print "The fetch program time: ", timeit.default_timer() - start_time

        # init article list
        article_lst.append(article_id)

    return(article_lst)

# matching articles
def matching(articles, keywords):
    weights = []
    combos = []
    combos_lst = []
    weighted_el_lst = []
    weighted_el_ = []
    for article in range(len(articles)):
        #if len(list(set(article) & set(article + 1))) != 0:
        if article is not len(articles) - 1 : # does not match the last index
            for index in range(len(articles)):
                if (index + article) <= len(articles) - 1: # article number + index does not equal the full length of articles
                    #print(article, index)
                    if len(list(set(articles[article]).intersection(set(articles[article + index])))) != 0: # if there are any matches
                        #empty list to make weighted edges form: (edges, edge, weights)
                        combos_lst = []
                        weights_lst = []

                        combo = tuple([keywords[article]] + [keywords[article + index]])
                        combos.append(combo)

                        #combos for weighted edges
                        combo_lst = keywords[article]
                        combos_lst.append(combo_lst)
                        combo_lst2 = keywords[article + index]
                        combos_lst.append(combo_lst2)

                        weight = len(list(set(articles[article]).intersection(set(articles[article + index]))))
                        weights.append(weight)

                        weighted_el = zip(combos, weights)
                        weighted_el_.append(weighted_el)

                        #for weighted edges
                        weights_lst.append(weight)
                        weighted_el_tup = tuple(combos_lst + weights_lst)
                        weighted_el_lst.append(weighted_el_tup)

    return(weights, combos, weighted_el_, weighted_el_lst)

# make a weighted graph of the keywords with the keyword of interest
def keyword_graph(weighted_el_lst):
    G = nx.Graph()

    # specify the keyword that you're interested in
    interested_keyword = raw_input("Please enter the keyword you are interested in: ")
    print "You enter: ", interested_keyword

    weighted_edgelist = match[3]

    # only choose the combos that include the interested keyword
    graph_of_interest = [interest for interest in weighted_edgelist if interested_keyword in interest]

    G.add_weighted_edges_from(graph_of_interest)
    pos = nx.spring_layout(G)
    custom_labels = {}
    nodes = G.nodes()
    for labels in nodes:
        custom_labels[labels] = labels

    edge_weight=dict([((u,v,),int(d['weight'])) for u,v,d in G.edges(data=True)])

    nx.draw_networkx_edge_labels(G,pos,edge_labels=edge_weight)
    nx.draw_networkx_nodes(G,pos)
    nx.draw_networkx_edges(G,pos)
    nx.draw_networkx_labels(G,pos)
    plt.show()


#############
# Extra code:

# Make an edgelist author (for a list of list)(Max)
def create_edgelist_author(author_goods, author_new):

    # creating a tuple
    author_combined = zip(author_goods, author_new)

    # edgelist code
    author_to_articles = {}
    for combo in author_combined:
        authors = combo[0]
        #author_list.append(authors)
        article_id = combo[1]
        #article_list.append(article_id)

        for author in authors:
            if author not in author_to_articles:
                author_to_articles[author] = [] #article_id goes in there
            #else:
            author_to_articles[author].append(article_id)

    return(author_to_articles)

# Make an edgelist author (for a list of list)(Max)
def create_edgelist_author(author_goods, author_new):

    # creating a tuple
    author_combined = zip(author_goods, author_new)

    # edgelist code
    author_to_articles = {}
    for combo in author_combined:
        authors = combo[0]
        #author_list.append(authors)
        article_id = combo[1]
        #article_list.append(article_id)

        for author in authors:
            if author not in author_to_articles:
                author_to_articles[author] = [] #article_id goes in there
            #else:
            author_to_articles[author].append(article_id)

    return(author_to_articles)

# Make a weighted graph (for a list of list) (max)
def create_weighted_graph(author_to_articles):
    G = nx.Graph()
    for author in author_to_articles:
        new_edges = list(itertools.combinations(author_to_articles[author], 2))
        for new_e in new_edges:
            if not G.has_edge(new_e[0], new_e[1]):
                G.add_edge(new_e[0], new_e[1], weight=1)
            else:
                G[new_e[0]][new_e[1]]["weight"] += 1
        G.edges(data=True)

    return(G)
