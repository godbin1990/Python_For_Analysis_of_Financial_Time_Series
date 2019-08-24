#coding=utf-8

import requests
import os
from bs4 import BeautifulSoup


base_url = r'https://faculty.chicagobooth.edu/ruey.tsay/teaching/fts3/'

def get_list():
    r=requests.get(base_url)
    soup=BeautifulSoup(r.content,'html5lib')
    links=soup.findAll('a')
    file_names=[link.string for link in links if link['href'][:4]=='file']
    for ii,filename in enumerate(file_names):
        start=0
        for index,i in enumerate(filename):
            if i!=' ':
                start=index
                break
        file_names[ii]=filename[start:]
    file_links=[base_url+name for name in file_names]
    return file_links

def request_url(strurl):
    proxy = '127.0.0.1:8087'

    local_proxies = {
        'http': 'http://' + proxy,
        'https': 'https://' + proxy,
    }

    try:
        r = requests.get(strurl)
    except:
        try:
            r = requests.get(strurl, proxies=local_proxies)
        except:
            r=None
    return r

def download():
    try:
        downloadlist=get_list()
    except:
        print('列表获取失败')
        return
    success=0
    print('下载开始,共%d个文件'%len(downloadlist))
    for i,file_url in enumerate(downloadlist):
        file_name=file_url.split('/')[-1]
        print('进度%d%%' % int(100 * i / len(downloadlist)))
        print('%s下载开始' % file_url)
        if not os.path.exists(file_name):
            try:
                for try_count in range(3):
                    rt=request_url(file_url)
                    if rt!=None:
                        break
                with open(file_name,'wb') as f:
                    f.write(rt.content)
                    f.close()
                    print('%s下载完成'%file_url)
                    success+=1
            except:
                print('%s下载失败' % file_url)
        else:
            success += 1
            print('%s已存在'%file_name)
        print('进度%d%%'%int(100*(i+1)/len(downloadlist)))
    print('下载结束，共%d个文件，成功%d个，失败%d个'%(len(downloadlist),success,len(downloadlist)-success))

download()