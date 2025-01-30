# Pubmed_Data_Analysis

Data binnentrekken

EUtilities
	ESearch => ID binnentrekken via zoekquery in URL mee te geven
	EFetch => Lijst van ID's gebruiken om data van alle artikelen binnen te trekken
	Met API key: 10requests/seconde, zonder 3requests/seconde
	Limiet op 10000 artikels om binnen te trekken => history server gebruiken
	History server => ID's kunnen worden opgeslagen als batch, gebruiken om data binnen te trekken.
	Volgende batch kan in history server worden opgeslagen enzovoort.History server "onthoudt" welke artikels al zijnbinnengehaald. Dus de volgende 10000 artikels zijn anderen

If you try to fetch records in batches without the History Server, you would have to manually run an ESearch query multiple times, each time requesting different slices of records.

For example, you might first run a search for articles 1-10,000, then for 10,001-20,000, and so on.

Problem: Each time you run ESearch, PubMed will return a fresh search result, but it doesn’t know where the last query left off. So, you’ll be re-running the same search over and over again, which is inefficient and will require more time and resources.

With the History Server, once you’ve done the search, you can continue fetching from where you left off by just referencing the WebEnv and QueryKey. You don’t have to repeat the search.


Edirect
	Via command line data binnentrekken
	Eenvoudiger maar minder flexibel
	Via docker container om te isoleren

