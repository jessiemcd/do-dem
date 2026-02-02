import numpy as np
import pickle
import pandas as pd
import astropy.time
import datetime
from datetime import timezone


def stdv_remove(stdv_file='/Users/jmdunca2/do-dem/reference_files/stdv_flares.pickle'):

    """
    Sometimes the stdv process identifies non-flares. Here we remove those from the list.
    
    """

    with open(stdv_file, 'rb') as f:
        data = pickle.load(f)

    firsts = [d[0] for d in data['stdv_flares']]
    #print(firsts)
    
    removing_these = [datetime.datetime(2014,12,11,19,3).replace(tzinfo=timezone.utc),
                      datetime.datetime(2014,12,11,19,6).replace(tzinfo=timezone.utc),
                      
                      datetime.datetime(2016,2,19,23,50).replace(tzinfo=timezone.utc),
                      
                      datetime.datetime(2017,10,10,4,45).replace(tzinfo=timezone.utc),
                      datetime.datetime(2017,10,10,5,11).replace(tzinfo=timezone.utc),
                      datetime.datetime(2017,10,10,5,13).replace(tzinfo=timezone.utc),
                      datetime.datetime(2017,10,10,5,16,50).replace(tzinfo=timezone.utc),
                      datetime.datetime(2017,10,10,5,18).replace(tzinfo=timezone.utc),
    
                      datetime.datetime(2018,5,29,21,35).replace(tzinfo=timezone.utc),

                      datetime.datetime(2018,9,7,15,45).replace(tzinfo=timezone.utc),
                      datetime.datetime(2018,9,7,15,48).replace(tzinfo=timezone.utc),
                      datetime.datetime(2018,9,7,15,50).replace(tzinfo=timezone.utc),
                      datetime.datetime(2018,9,7,15,51).replace(tzinfo=timezone.utc),
                      datetime.datetime(2018,9,7,19,6).replace(tzinfo=timezone.utc),
                      
                      datetime.datetime(2020,1,29,16,50).replace(tzinfo=timezone.utc),
                      datetime.datetime(2020,1,29,17,53,30).replace(tzinfo=timezone.utc),
                      
                      datetime.datetime(2021,1,8,10,14,30).replace(tzinfo=timezone.utc),
                      datetime.datetime(2021,1,14,13,27,10).replace(tzinfo=timezone.utc),
                      datetime.datetime(2021,11,17,18,23).replace(tzinfo=timezone.utc),
                      datetime.datetime(2021,11,19,23,15).replace(tzinfo=timezone.utc),
                      datetime.datetime(2021,11,22,0,42,35).replace(tzinfo=timezone.utc),
                      datetime.datetime(2022,6,3,22,23).replace(tzinfo=timezone.utc)
                      
                     ]

    flares = data['stdv_flares']
    removals = []
    for i in range(0, len(flares)):
        for rt in removing_these:
            if np.logical_and(rt > flares[i][0], rt < flares[i][1]):
                removals.append(i)
                print('removing: ', flares[i][0].strftime('%D %H-%M-%S'), flares[i][1].strftime('%D %H-%M-%S'))
                
    
    keepflares = [flares[i] for i in range(0, len(flares)) if i not in removals]
    data['stdv_flares'] = keepflares
    
    with open(stdv_file, 'wb') as f:
         # Pickle the 'data' dictionary using the highest protocol available.
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL) 

    
        
    


def show_manual_flares(manual_file='/Users/jmdunca2/do-dem/reference_files/manual_flares.pickle'):

    """
    Makes a giant list of manual flare times + adds them to the relevant file, prints a list. 
    See function itself to add a new manual flare. 
    """
    
    mfts = []
    mfts.append([astropy.time.Time('2016-04-22 20:40:00'), astropy.time.Time('2016-04-22 20:45:00')])
    mfts.append([astropy.time.Time('2016-04-22 22:35:00'), astropy.time.Time('2016-04-22 22:45:00')])
    
    mfts.append([astropy.time.Time('2016-07-27 00:29:30'), astropy.time.Time('2016-07-27 01:00:00')])
    
    mfts.append([astropy.time.Time('2017-09-11 15:55:00'), astropy.time.Time('2017-09-11 17:00:00')])
    mfts.append([astropy.time.Time('2017-09-11 17:30:00'), astropy.time.Time('2017-09-11 18:40:00')])
    mfts.append([astropy.time.Time('2017-09-11 19:00:00'), astropy.time.Time('2017-09-11 20:20:00')])
    mfts.append([astropy.time.Time('2017-09-11 20:40:00'), astropy.time.Time('2017-09-11 20:48:00')])
    
    mfts.append([astropy.time.Time('2017-09-12 22:58:00'), astropy.time.Time('2017-09-12 23:15:00')])
    mfts.append([astropy.time.Time('2017-09-12 21:42:00'), astropy.time.Time('2017-09-12 21:50:00')])
    
    
    mfts.append([astropy.time.Time('2018-05-29 16:01:00'), astropy.time.Time('2018-05-29 16:03:00')])
    mfts.append([astropy.time.Time('2018-05-29 16:10:00'), astropy.time.Time('2018-05-29 16:14:00')])
    mfts.append([astropy.time.Time('2018-05-29 17:34:00'), astropy.time.Time('2018-05-29 17:40:00')])
    mfts.append([astropy.time.Time('2018-05-29 19:19:00'), astropy.time.Time('2018-05-29 19:25:00')])
    
    mfts.append([astropy.time.Time('2018-05-29 21:31:45'), astropy.time.Time('2018-05-29 21:38:00')])
    
    mfts.append([astropy.time.Time('2019-04-12 17:12:00'), astropy.time.Time('2019-04-12 17:25:00')])
    mfts.append([astropy.time.Time('2019-04-12 18:20:00'), astropy.time.Time('2019-04-12 18:35:00')])
    
    mfts.append([astropy.time.Time('2019-04-13 04:21:00'), astropy.time.Time('2019-04-13 04:40:00')])
    mfts.append([astropy.time.Time('2019-04-13 06:10:00'), astropy.time.Time('2019-04-13 06:50:00')])
    mfts.append([astropy.time.Time('2019-04-13 08:45:00'), astropy.time.Time('2019-04-13 10:10:00')])
    
    mfts.append([astropy.time.Time('2020-06-06 19:35:00'), astropy.time.Time('2020-06-06 19:45:00')])
    mfts.append([astropy.time.Time('2020-06-06 19:55:00'), astropy.time.Time('2020-06-06 20:10:00')])
    
    mfts.append([astropy.time.Time('2020-06-07 19:49:00'), astropy.time.Time('2020-06-07 19:58:00')])
    mfts.append([astropy.time.Time('2020-06-07 21:25:30'), astropy.time.Time('2020-06-07 21:30:00')])
    
    mfts.append([astropy.time.Time('2020-06-08 21:00:00'), astropy.time.Time('2020-06-08 21:05:00')])
    mfts.append([astropy.time.Time('2020-06-08 21:50:00'), astropy.time.Time('2020-06-08 22:05:00')])
    
    mfts.append([astropy.time.Time('2021-04-29 21:25:00'), astropy.time.Time('2021-04-29 21:40:00')])
    
    mfts.append([astropy.time.Time('2021-05-03 17:00:00'), astropy.time.Time('2021-05-03 17:08:30')])
    mfts.append([astropy.time.Time('2021-05-03 17:37:00'), astropy.time.Time('2021-05-03 17:50:00')])
    
    mfts.append([astropy.time.Time('2021-07-20 08:47:00'), astropy.time.Time('2021-07-20 09:00:00')])
    mfts.append([astropy.time.Time('2021-07-20 10:00:00'), astropy.time.Time('2021-07-20 10:07:00')])
    mfts.append([astropy.time.Time('2021-07-20 10:57:35'), astropy.time.Time('2021-07-20 10:57:59')])
    mfts.append([astropy.time.Time('2021-07-20 11:45:00'), astropy.time.Time('2021-07-20 11:50:00')])
    
    mfts.append([astropy.time.Time('2021-11-20 00:25:00'), astropy.time.Time('2021-11-20 00:40:00')])
    
    mfts.append([astropy.time.Time('2022-09-06 17:20:00'), astropy.time.Time('2022-09-06 18:10:00')])
    mfts.append([astropy.time.Time('2022-09-06 18:50:00'), astropy.time.Time('2022-09-06 19:50:00')])
    
    mfts.append([astropy.time.Time('2022-12-09 23:39:00'), astropy.time.Time('2022-12-09 23:46:00')])
    mfts.append([astropy.time.Time('2022-12-09 23:47:00'), astropy.time.Time('2022-12-09 23:55:00')])
    mfts.append([astropy.time.Time('2022-12-09 23:55:00'), astropy.time.Time('2022-12-10 00:00:30')])
    mfts.append([astropy.time.Time('2022-12-10 00:30:00'), astropy.time.Time('2022-12-10 00:45:00')])

    print('')
    print('Here are the manual flares; edit the function if you want more.')
    for m in mfts:
        print(m[0].strftime('%D %H-%M-%S'), m[1].strftime('%D %H-%M-%S'))
    
    
    data = {'manual_flares': mfts}
    
    with open(manual_file, 'wb') as f:
         # Pickle the 'data' dictionary using the highest protocol available.
         pickle.dump(data, f, pickle.HIGHEST_PROTOCOL) 