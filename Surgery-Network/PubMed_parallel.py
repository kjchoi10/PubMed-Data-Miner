### testing parallel computing
import multiprocessing as mp
import random
import string
from multiprocessing.dummy import Pool as ThreadPool
from Bio import Entrez
from Bio import Medline

## pooling example 1/21/16
from multiprocessing.dummy import Pool as ThreadPool
from Bio import Entrez
from Bio import Medline

class setup():
    # setup the email
    def email(self):
        email = raw_input("Please enter your email: ")
        print "You enter: ", email
        return(email)

init = setup()
Entrez.email = init.email()

pool = ThreadPool(4)

item = ['9749847','9755759','9757885','9766300','9769004','9784814','9800512','9804224','9877391','9879865','9888139','9918973','9919913']

record = Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs",from_uid=item)

handle = Entrez.efetch(db="pubmed", id=item, rettype="medline",retmode="xml")

results = pool.map(Entrez.read(handle), handle)

pool.close()
pool.join()











### parallel using pool

# list of pmid ids
lst = ['9749847','9755759','9757885','9766300','9769004'],'9784814','9800512','9804224','9877391','9879865','9888139','9918973','9919913']
lst = ['9749847']

lst =  [['9749847','9755759'],['9757885','9766300','9769004'],['9784814','9800512','9804224','9877391'],['9879865','9888139','9918973','9919913']]

lst = ['9749847','9755759','9757885','9766300','9769004','9784814','9800512','9804224','9877391','9879865','9888139','9918973','9919913']

# define function
def list_records(item):
    record = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs",from_uid=item))
    #print(item)
    return(record)

pool = mp.Pool(processes=4)
results = [pool.apply(list_records, args=(lst, ))]
print(results)
output = [p.get() for p in results]
print(output)


# define a example function
def rand_string(length, output):
    """ Generates a random string of numbers, lower- and uppercase chars. """
    rand_str = ''.join(random.choice(
                    string.ascii_lowercase
                    + string.ascii_uppercase
                    + string.digits)
               for i in range(length))
    output.put(rand_str)

# Setup a list of processes that we want to run
processes = [mp.Process(target=rand_string, args=(5, output)) for x in range(4)]

# Run processes
for p in processes:
    p.start()

# Exit the completed processes
for p in processes:
    p.join()

# Get process results from the output queue
results = [output.get() for p in processes]

print(results


# list of records from the splitted pmid ids: this is the right one -- use parallel computing for this
new = split_pmid[0]
def list_records(split_pmid):
    record_list = []
    for splits in range(len(split_pmid)):
        list_pmid = split_pmid[splits]
        temp = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs",from_uid=list_pmid))
        record_list.append(temp)

    return(record_list)

# try threading instead of paralleing computing
import threading
import Queue
import time

# input queue to be processed by many threads
q_in = Queue.Queue(maxsize=0)

# output queue to be processed by one thread
q_out = Queue.Queue(maxsize=0)

# number of worker threads to complete the processing
num_worker_threads = 20

# the lst
lst =  ['9749847','9755759','9757885','9766300','9769004','9784814','9800512','9804224','9877391','9879865','9888139','9918973','9919913']

# queues data strucutre only accepts elements and not lists

# another test code:
# try threading instead of paralleing computing
import threading
import Queue
import time

# input queue to be processed by many threads
q_in = Queue.Queue(maxsize=0)

# output queue to be processed by one thread
q_out = Queue.Queue(maxsize=0)

# number of worker threads to complete the processing
num_worker_threads = 20

# the lst
lst =  ['9749847','9755759','9757885','9766300','9769004','9784814','9800512','9804224','9877391','9879865','9888139','9918973','9919913']

# queues data strucutre only accepts elements and not lists

# another test code:
def worker():
    item =


# process that each worker thread will execute until the Queue is empty
def worker():
    while True:
        # get item from queue, do work on it, let queue know processing is done for one item
        item = q_in.get()
        q_out.put(list_records(item))
        q_in.task_done()

def list_records(item): # pass in a list len 100
    results = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs",from_uid=item))
    return(results)

# another queued thread we will use to print output
def printer():
    while True:
        # get an item processed by worker threads and print the result. Let queue know item has been processed
        item = q_out.get()
        print("step")
        q_out.task_done()

# launch all of our queued processes
def main():
    # Launches a number of worker threads to perform operations using the queue of inputs
    for i in range(num_worker_threads):
         t = threading.Thread(target=worker)
         #t.daemon = True
         threads.append(t)
         t.start()

    # launches a single "printer" thread to output the result (makes things neater)
    #t = threading.Thread(target=printer)
    #t.daemon = True
    #t.start()

    # put items on the input queue (numbers to be squared)
    for item in lst:
        q_in.put(item)

    # wait for two queues to be emptied (and workers to close)
    q_in.join()       # block until all tasks are done
    q_out.join()

    print "Processing Complete"


# process that each worker thread will execute until the Queue is empty
def worker():
    while True:
        # get item from queue, do work on it, let queue know processing is don  e for one item
        item = q_in.get()
        q_out.put(list_records(item))
        q_in.task_done()

def list_records(item): # pass in a list len 100
    results = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs",from_uid=item))

    return(results)

# another queued thread we will use to print output
def printer():
    while True:
        # get an item processed by worker threads and print the result. Let queue know item has been processed
        item = q_out.get()
        print("step")
        q_out.task_done()

# launch all of our queued processes
def main():
    # Launches a number of worker threads to perform operations using the queue of inputs
    for i in range(num_worker_threads):
         t = threading.Thread(target=worker)
         #t.daemon = True
         threads.append(t)
         t.start()

    # launches a single "printer" thread to output the result (makes things neater)
    #t = threading.Thread(target=printer)
    #t.daemon = True
    #t.start()

    # put items on the input queue (numbers to be squared)
    for item in lst:
        q_in.put(item)

    # wait for two queues to be emptied (and workers to close)
    q_in.join()       # block until all tasks are done
    q_out.join()

    print "Processing Complete"

main()



def list_records(worker, split_pmid):
    record_list = []
    for splits in range(len(split_pmid)):
        list_pmid = split_pmid[splits]
        temp = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs",from_uid=list_pmid))
        record_list.append(temp)

    return(record_list)

def threader():
    while True:
        split_pmid = q.get()
        list_records(worker, split_pmid)
        q.task_done()

q = Queue()

for x in range(len(split_pmid)): # about 100
    t = threading.Thread(target = threader)
    t.daemon = True
    t.start()

start = time.time()

for worker in range(100):
    q.put(worker)

q.join()

print("time took: ", time.time() - start)

##### end of test

### code from batching in PubMed_good.py


        # logic behind extracting data from medline
        if len(pmid) <= 10000:
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline",retmode="text")
            parse = Medline.parse(handle)
            records = list(parse)

        number = len(pmid)/10000
        print number
        # splits the pmid list into a size that efetch can work with (10,000 upper limit)
        split_pmid = list(split_every(number, pmid))
        return(records, split_pmid)

        else:
            # find divisor that will make len(pmid) ~ 10,000
            number = len(pmid)/10000
            print number
            # splits the pmid list into a size that efetch can work with (10,000 upper limit)
            split_pmid = list(split_every(number, pmid))

        return(split_pmid)

            records = []
            # extract data
            for item in range(len(split_pmid)): # trying doing range instead of in
                handle = Entrez.efetch(db="pubmed", id=item, rettype="medline",retmode="text")
                parse = Medline.parse(handle)
                records.append(parse)


### batch code:
# Getting data from Medline: journals, authors, location
def record_fetch(pmid, webenv, query_key):
    # length of pmid - test sake
    print "The length of pmid is", len(pmid)

    # timer

    # use session history to batch queries
    # timer
    start_time_1 = timeit.default_timer()
    start = 0
    batch_size = len(pmid)
    fetch_handle = Entrez.efetch(db="pubmed", rettype="medline", retmode="text", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key)
    print "The handle fetch program time is ", timeit.default_timer() - start_time_1

    start_time_1 = timeit.default_timer()
    count = len(pmid)
    batch_size = 1000
    for start in range(0, count, batch_size):
        end = min(count, start+batch_size)
        print("Going to download record %i to %i" % (start+1, end))
        fetch_handle = Entrez.efetch(db="pubmed", rettype="medline", retmode="text", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key)
        data = fetch_handle.read()
        fetch_handle.close()
    print "The handle fetch program time is ", timeit.default_timer() - start_time_1

    ## with attemp and except
    start_time_1 = timeit.default_timer()
    count = len(pmid)
    batch_size = 1000
    try:
        from urllib.error import HTTPError # for Python 3
    except ImportError:
        from urllib2 import HTTPError # for Python 2
    for start in range(0, count, batch_size):
        end = min(count, start+batch_size)
        print("Going to download record %i to %i" % (start+1, end))
        attempt = 1
        while attempt <= 3:
            try:
                fetch_handle = Entrez.efetch(db="pubmed", rettype="medline", retmode="text", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key)
            except HTTPError as err:
                if 500 <= err.code <= 599:
                    print("Received error from server %s" % err)
                    print("Attempt %i of 3" % attempt)
                    attempt += 1
                    time.sleep(2)
                else:
                    raise
        records = Medline.parse(fetch_handle)
        records = list(records)
        fetch_handle.close()
    print "The handle fetch program time is ", timeit.default_timer() - start_time_1

    # use Medline parser to turn data into a parser object
    start_time_2 = timeit.default_timer()
    records = Medline.parse(fetch_handle)
    records = list(records)
    print "The Medline parse program time is ", timeit.default_timer() - start_time_2


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
    print "The author program time is ", timeit.default_timer() - start_time_3

    # data frame
    data = pd.DataFrame({'authors': record_authors, 'firstauthor': record_first, 'journals': record_journals, 'datepublication': record_dp,
    'place': record_place, 'mesh': record_mesh, 'pmid': record_pmid,'affliation': record_aff}, columns=['authors', 'firstauthor', 'journals',
    'datepublication', 'place', 'mesh', 'pmid', 'affliation'])

    return(data)


        # logic behind extracting data from medline
        if len(pmid) <= 10000:
            handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline",retmode="text")
            parse = Medline.parse(handle)
            records = list(parse)

        number = len(pmid)/10000
        print number
        # splits the pmid list into a size that efetch can work with (10,000 upper limit)
        split_pmid = list(split_every(number, pmid))
        return(records, split_pmid)

        else:
            # find divisor that will make len(pmid) ~ 10,000
            number = len(pmid)/10000
            print number
            # splits the pmid list into a size that efetch can work with (10,000 upper limit)
            split_pmid = list(split_every(number, pmid))

        return(split_pmid)

            records = []
            # extract data
            for item in range(len(split_pmid)): # trying doing range instead of in
                handle = Entrez.efetch(db="pubmed", id=item, rettype="medline",retmode="text")
                parse = Medline.parse(handle)
                records.append(parse)
