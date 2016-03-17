### Libraries
from Bio import Entrez
from Bio import Medline
from collections import Counter
from collections import defaultdict
import collections
import itertools
import timeit
import time
from itertools import islice
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from mpl_toolkits.basemap import Basemap
from geopy.geocoders import Nominatim

# essential Libraries
from Bio import Entrez
from Bio import Medline
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import timeit
from collections import defaultdict


###

# Execution code

### 1: Setup

# Best version of code so far (automatic)
init = setup()
Entrez.email = init.email()
init_keyword = init.keyword()
keywords = ["Robotic Surgery", "Cardiac Surgery", "Colorectal Surgery", "General Surgery", "Gynecologic Surgery", "Head & Neck Surgery","Thoracic Surgery", "Urologic Surgery"]
keywords2 = ['Adolescent Medicine', 'Aerospace Medicine', 'Allergy and Immunology', 'Andrology', 'Anesthesiology', 'Bariatric Medicine', 'Behavioral Medicine', 'Clinical Medicine', 'Community Medicine', 'Dermatology',
            'Disaster Medicine', 'Emergency Medicine', 'Forensic Medicine', 'General Practice', 'Genetics, Medicine', 'Geography, Medicine', 'Geriatrics, Medicine', 'Global Health', 'Hospital Medicine', 'Integrative Medicine',
            'Internal Medicine', 'Military Medicine', 'Molecular Medicine', 'Naval Medicine', 'Neurology', 'Osteopathic Medicine', 'Palliative Medicine', 'Pathology', 'Pediatrics', 'Physical and Rehabilitation Medicine',
            'Psychiatry', 'Public Health', 'Radiology', 'Regenerative Medicine', 'Reproductive Medicine', 'Social Medicine', 'Specialties, Surgical', 'Sports Medicine', 'Telemedicine', 'Theranostic Nanomedicine',
            'Travel Medicine', 'Tropical Medicine', 'Venerology', 'Wilderness Medicine', 'Mortuary Practice', 'Nursing', 'Nursing, Practical', 'Nutritional Sciences', 'Optometry', 'Orthoptics', 'Pharmacology', 'Pharmacy',
            'Podiatry', 'Psychology, Medical', 'Serology', 'Sociology, Medical', 'Toxicology', 'Veterinary Medicine']
keywords3 = ['Urologic Surgical Procedures','Cystectomy','Cystoscopy','Cystotomy','Kidney Transplantation','Nephrectomy','Ureteroscopy','UrinaryDiversion','Cystostomy','Nephrostomy, Percutaneous','Ureterostomy','Urologic Surgical Procedures, Male','Circumcision, Male','Orchiectomy','Orchiopexy','Penile Implantation','Prostatectomy','Vasectomy','Vasovasostomy']
df = analysis(keywords, date)
articles = [item['pmid'] for item in df]
place = [item['place'] for item in df]
date = [item['datepublication'] for item in df]
affliation = [item['affliation'] for item in df]
match_articles = matching(articles, keywords)
match_place = matching(place, keywords)
match_date = matching(date, keywords)
match_affliation = matching(affliation, keywords)
graph = target_graph(match[3])
graph = all_graph(match[3])

#### plotting graphs:
# Place for first keyword
place = df[0]['place']
place = map(lambda x: x.upper(), place)
df[0]['place'] = map(lambda x: x.upper(), df[0]['place']) # to avoid duplicates of country names
df[0]['affliation'] = map(lambda x: x.upper(), df[0]['affliation'])
df['journals'] = map(lambda x: x.upper(), df['journals'])
place = df[0]['place']
affliation = df[0]['affliation']
journal = df['journals']
p1 = place
p2 = affliation
p3 = journal
vc = p1.value_counts()
vc2 = p2.value_counts()
vc3 = p3.value_counts()
vc = vc.sort_index()
vc2 = vc2.sort_index
vc3 = vc3.sort_index
vc.plot(kind='bar')
unique_countries = pd.unique(df[0].place.ravel()) # getting the unique country names
unique_affliations = pd.unique(df[0].affliation.ravel()) # getting the unique country names
unique_journals = pd.unique(df.journal.ravel())

# for urologic surgeryp1 = place[0]
p1 = place[6]
vc = p1.value_counts()
vc = vc.sort_index()
vc.plot(kind='bar')
