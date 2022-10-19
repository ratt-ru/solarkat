from os import sys
from astropy.time import Time
import numpy as numpy

    
def split_times(interval,ms):
  dic={ "-": "/", " ":"/"}

  def replace_all(text, dic):
      
      for i, j in dic.items():
          text = text.replace(i, j)
      return text

  array=[]
  tt = tb.open(ms)
  all_times = list(numpy.unique(tb.getcol('TIME')))
  t0 = all_times[0]
  t1 = all_times[-1]
  dt = (t1-t0)/(interval)
  for i in range(interval):
      t2=dt*i+t0
      t_iso = Time(t2/86400.0,format='mjd').iso
      array.append(t_iso)
  
  for i in range(len(array)):
      SPLITED_ms=ms.replace('.ms','_time'+str(i)+'.ms')
      split(vis=ms,outputvis='/vault-ike/kincaid/solarkat'+SPLITED_ms, timerange=replace_all(array[i],dic)+'~'+replace_all(array[i+1],dic) ,datacolumn='all')
      print(SPLITED_ms + 'Done')


def split_scans(ms):

  maintab = tb.open(ms)
  scans = list(numpy.unique(tb.getcol('SCAN_NUMBER')))

  for scan in scans:
      SPLITED_ms=ms.replace('.ms','_scan_'+str(scan)+'.ms')
      split(vis=ms,outputvis='/vault-ike/kincaid/solarkat/splitted_ms/'+SPLITED_ms, scan=str(scan) ,datacolumn='all')
      print(SPLITED_ms + 'Done')


if __name__ == "__main__":
       ms=sys.argv[3]
       interval=30
       #split_times(interval,ms)
       split_scans(ms)
#

  
