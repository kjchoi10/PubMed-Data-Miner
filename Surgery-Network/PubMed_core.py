# matching, analysis, and graph algorithms

# run analysis and create a list of suitable article ids - this code takes way too long - why? much longer than when i execute indiv
# number of days in the past you want to look at
date = 1825
# 5 years = 1825
# 4 years = 1460
# 3 years = 1095
# 2 years = 730
# 1 year = 365

# min and max date range:
date = ['2015/01/01', '2014/01/01', '2013/01/01', '2012/01/01', '2011/01/01', '2010/01/01']

df_0.7 = df.sample(frac=0.7)


# matching, analysis, and graph algorithms
def analysis(keywords, date):
    article_lst = []
    ############################
    # [ initial handler setup ]
    start_time_final = timeit.default_timer()
    for word in range(len(keywords)):
        print(keywords[word])
        # to find keyword population with count
        start_time = timeit.default_timer()
        #search_handle = Entrez.esearch(db="pubmed",term=keywords[word],retmax=1, reldata=date, usehistory="y")
        search_handle = Entrez.esearch(db="pubmed",term=keywords[word],retmax=1, mindate='2015', maxdate='2014', usehistory="y")
        #handle = Entrez.esearch(db="pubmed",term=keywords[0], retmax=1)
        #handle = Entrez.esearch(db="pubmed",term=keywords[1],retmax=1)
        search_results = Entrez.read(search_handle)
        count = search_results['Count']

        if int(count) <= 10000:
            sample = 10000
        else:
            sample = float(count) * 0.30
            sample = int(sample)

        ### 1. I want to download the full sample of data (will take a long time) for the specified year
        # 2. Then randomly select 1000 articles for that year
        # 3. need to get the min and max date more automatic than changing manually

        # 1. Date specified

        # taking 10% of population and using this as sample
        #search_handle= Entrez.esearch(db="pubmed",term=keywords[word],retmax=sample, reldata=date, usehistory="y")
        search_handle = Entrez.esearch(db="pubmed",term=keywords[word],retmax=sample, mindate='2015', maxdate='2014', usehistory="y")
        search_results = Entrez.read(search_handle)

        # web servery history parameters
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]

        print ("The initial handle program time: ", timeit.default_timer() - start_time)

        ############################
        # [ handle test data ]
        #handle = Entrez.esearch(db="pubmed",term=keywords[0],retmax=sample)
        #handle = Entrez.esearch(db="pubmed",term=keywords[1],retmax=sample)

        print ("Sample size is ", sample)

        # storing the sample as the idlist for the articles

        ############################
        # [ splitting step ]
        start_time_2 = timeit.default_timer()
        id_results = search_results["IdList"]
        pmid = id_results
        pmid = sorted(pmid)
        print ("The sort program time ", timeit.default_timer() - start_time_2)
        test = record_fetch(pmid)
        article_lst.append(test)
    print ("The final program time: ", timeit.default_timer() - start_time_final)
    return(article_lst)

# matching articles
def matching(df, keywords):
    articles = [item['pmid'] for item in df]
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
def target_graph(weighted_el_lst):
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

# make a weighted graph of the keywords
def all_graph(weighted_el_lst):
    G = nx.Graph()
    preds = nx.adamic_adar_index(G, match[1])

    weighted_edgelist = match[3]
    G.add_weighted_edges_from(weighted_edgelist)
    preds = nx.adamic_adar_index(G, match[1])
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
