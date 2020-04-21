#! python
# -*- coding: utf-8 -*-
# author : charlie

import os
import re
import time
import requests
from playsound import playsound
from bs4 import BeautifulSoup as bs


word = input('$ -> ').strip()
while word != 'exit':
	if word == "":
		word = input('$ -> ').strip()
		continue
	url = 'https://cn.bing.com/dict/search?q=' + word + '&FORM=HDRSC6'
	headers = {
		'cookie':'_EDGE_S=F=1&SID=32BCC069762562CB21D8CD28770B631F; _EDGE_V=1; MUID=263CE2A7B327662C3E05EFE6B209670A; MUIDB=263CE2A7B327662C3E05EFE6B209670A; SNRHOP=I=&TS=; SRCHD=AF=NOFORM; SRCHUID=V=2&GUID=74052B20B3E94D2785530401B7A3903D&dmnchg=1; SRCHUSR=DOB=20190417&T=1555512247000; SRCHHPGUSR=CW=1582&CH=822&DPR=1.2000000476837158&UTC=480&WTS=63691109047; _SS=SID=32BCC069762562CB21D8CD28770B631F&HV=1555512254',
		'referer':'https://cn.bing.com/dict?FORM=HDRSC6',
		'upgrade-insecure-requests':'1',
		'user-agent':'Mozilla/5.0 (Windows NT 10.0; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/63.0.3239.132 Safari/537.36',
	}
	resp = requests.get(url, headers = headers)
	soup = bs(resp.text, 'lxml')
	try:
		# 输出含义
		hd = soup.find('div', class_ = 'hd_prUS').text + '  ' +  soup.find('div', class_ = 'hd_pr').text
		print(hd)
		defines = soup.find('div',class_ = 'qdef').find_all('li')
		for define in defines:
			pos = define.find('span', class_ = 'pos')
			def_ = define.find('span', class_ = 'def')
			if '网络' in pos.text:
				print('network：' + def_.text)
			else:
				print(pos.text + ' ' + def_.text)

		# 输出例句
		exams = soup.find('div', class_ = 'de_seg').find_all('div', class_ = 'se_lis')
		print('eg.')
		for exam in exams:
			print(exam.text)

		# 语音
		try: 
			voice = soup.find('a', class_ = 'bigaud')['onclick']
			#print(voice)
			partten = re.compile(r'https?://[-A-Za-z0-9+&@#/%?=~_|!:,.;]+[-A-Za-z0-9+&@#/%=~_|]')
			voice_url = re.findall(partten, voice)
			# print(voice_url[0] )
			voice_response = requests.get( voice_url[0] )
			if voice_response.status_code == 200:
			    with open("bing_voice.mp3" , "wb") as code: 
			        code.write(voice_response.content)
			    playsound('bing_voice.mp3')
			    # print('removing ')
			    os.remove('bing_voice.mp3')
		except:
			print('sorry, no pronunciation available ~ ')

	except:
		# 处理中文
		examples = []
		try:
			defines = soup.find('div',class_ = 'qdef').find_all('li')
			for define in defines:
				pos = define.find('span', class_ = 'pos')
				def_ = define.find('span', class_ = 'def')
				if '网络' in pos.text:
					print('network：' + def_.text)
				else:
					print(pos.text + ' ' + def_.text)

			# 输出例句
			exams = soup.find('div', class_ = 'df_div').find_all('div', class_ = 'def_fl')
			for exam in exams:
				examples += exam.find_all('div', class_ = 'df_cr_w')
			print('汉英：', end = ' ')
			for example in examples:
				print(example.text)
			print()
		except:
			try:
				print('translation -> ', end = '')
				tran = soup.find('div', class_ = 'lf_area').find('div', class_ = 'p1-11')
				print(tran.text)
			except:
				print("     >……<   ")

	word = input('\n$ -> ').strip()

#os.system('exit')
