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
    print(array)
    
    if i < (len(array)-1):
      print(i)
      print(i+1)
      print(array[i])
      SPLITED_ms=ms.replace('.ms','_time'+str(i)+'.ms')
      split(vis=ms,outputvis='/vault-ike/kincaid/solarkat/'+SPLITED_ms, timerange=replace_all(array[i],dic)+'~'+replace_all(array[i+1],dic) ,datacolumn='all')
      print(str(SPLITED_ms), 'Done')
    else:
      print('Done')
      exit()


if __name__ == "__main__":
  ms=sys.argv[3]
  interval=int(sys.argv[4])
  split_times(interval,ms)
