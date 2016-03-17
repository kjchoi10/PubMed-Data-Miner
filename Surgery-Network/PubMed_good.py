# setting up with your email
class setup():
    # setup the email
    def email(self):
        email = raw_input("Please enter your email: ")
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


# Getting data from Medline: journals, authors, location
def record_fetch(pmid):
    # length of pmid - test sake
    print ("The length of pmid is", len(pmid))

    # timer
    start_time_1 = timeit.default_timer()

    # fetching pubmed articles using pmid ids
    fetch_handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline",retmode="text")

    # timer
    print ("The handle fetch program time is ", timeit.default_timer() - start_time_1)

    # timer
    start_time_2 = timeit.default_timer()

    # parsing Medline
    records = Medline.parse(fetch_handle)
    records = list(records)

    # timer
    print ("The Medline parse program time is ", timeit.default_timer() - start_time_2)

    # timer
    start_time_3 = timeit.default_timer()

    # init data columns
    record_authors = []
    record_first = []
    record_journals = []
    record_dp = []
    record_place = []
    record_mesh = []
    record_pmid = []
    record_aff = []

    # iterate over records - try and errors are there to catch None types
    for record in records:
        try:
            record_authors.append(record.get("AU"))
        except TypeError:
            record_authors.append("None")
        try:
            record_first.append(record.get("AU")[0])
        except TypeError:
            record_first.append("None")
        try:
            record_journals.append(record.get("JT"))
        except TypeError:
            record_journals.append("None")
        try:
            record_dp.append(record.get("DP"))
        except TypeError:
            record_dp.append("None")
        try:
            record_place.append(record.get("PL"))
        except TypeError:
            record_place.append("None")
        try:
            record_mesh.append(record.get("MH"))
        except TypeError:
            record_mesh.append("None")
        try:
            record_pmid.append(record.get("PMID"))
        except TypeError:
            record_pmid.append("None")
        try:
            record_aff.append(record.get("AD"))
        except TypeError:
            record_aff.append("None")

    # timer
    print ("The author program time is ", timeit.default_timer() - start_time_3)

    # data frame
    data = pd.DataFrame({'authors': record_authors, 'firstauthor': record_first, 'journals': record_journals, 'datepublication': record_dp,
    'place': record_place, 'mesh': record_mesh, 'pmid': record_pmid,'affliation': record_aff}, columns=['authors', 'firstauthor', 'journals',
    'datepublication', 'place', 'mesh', 'pmid', 'affliation'])

    # dropping miss data (None types)
    data = data.dropna()

    return(data)

# mapping  countries
def map_countries(data):
    geolocator = Nominatim()
    location = [geolocator.geocode(country) for country in unique_countries]
    #institution = [geolocator.geocode(affliation) for affliation in unique_affliations]
    print(location)
    # make sure the value of resolution is a lowercase L,
    #  for 'low', not a numeral 1
    my_map = Basemap(projection='robin', lat_0=0, lon_0=-100, resolution='l', area_thresh=1000.0)

    my_map.drawcoastlines()
    my_map.drawcountries()
    my_map.fillcontinents(color='coral')
    my_map.drawmapboundary()

    my_map.drawmeridians(np.arange(0, 360, 30))
    my_map.drawparallels(np.arange(-90, 90, 30))

    lons = [(location[lon][1][1]) for lon in range(len(location))]
    lats = [(location[lat][1][0]) for lat in range(len(location))]
    x,y = my_map(lons, lats)
    my_map.plot(x, y, 'bo', markersize=10)
    plt.show()
