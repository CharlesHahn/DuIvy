# author: charlie

#####################
# 1. show Greek alphabet done
# 2. amino searching 
#
#####################

from PIL import Image

amino_data = [
	"G ; Gly ; Glycine ; 甘氨酸",             "A ; Ala ; Alanine ; 丙氨酸",
	"V ; Val ; Valine ; 缬氨酸",              "L ; Leu ; Leucine ; 亮氨酸",
	"I ; Ile ; Isoleucine ; 异亮氨酸",        "S ; Ser ; Serine ; 丝氨酸", 
	"T ; Thr ; Threonine ; 苏氨酸",           "F ; Phe ; Phenylalanine ; 苯丙氨酸",
	"Y ; Tyr ; Tyrosine ; 络氨酸",            "W ; Trp ; Tryptophan ; 色氨酸",
	"D ; Asp ; Aspartate ; 天冬氨酸",         "E ; Glu ; Glutamate ; 谷氨酸",
	"N ; Asn ; Asparagine ; 天冬酰胺",        "Q ; Gln ; Glutamine ; 谷酰胺",
	"K ; Lys ; Lysine ; 赖氨酸",              "R ; Arg ; Arginine ; 精氨酸",
	"H ; His ; Histidine ; 组氨酸",           "C ; Cys ; Cysteine ; 半胱氨酸",
	"M ; Met ; Methionine ; 甲硫氨酸",        "P ; Pro ; Proline ; 脯氨酸"
	]

amino_infomation = """{:<10}{:<10}\n{:<10}{:<10}\n{:<10}{:<10}\n{:<10}{:<10}\n{:<10}{:<10}\n{:<10}{:<10}
""".format(	"脂肪", "G A V L I", "羟基", "S T", "苯环", "F Y W ", "酸性", "D E", 
	"碱性", "K R H ", "酰胺", "N Q", "硫", "C M ", "亚氨基酸", "P",)


greek = """{:<3}{:<3}{:<10}  {:<3}{:<3}{:<10}  {:<3}{:<3}{:<10}  {:<3}{:<3}{:<10}
{:<3}{:<3}{:<10}  {:<3}{:<3}{:<10}  {:<3}{:<3}{:<10}  {:<3}{:<3}{:<10}
{:<3}{:<3}{:<10}  {:<3}{:<3}{:<10}  {:<3}{:<3}{:<10}  {:<3}{:<3}{:<10}
{:<3}{:<3}{:<10}  {:<3}{:<3}{:<10}  {:<3}{:<3}{:<10}  {:<3}{:<3}{:<10}
{:<3}{:<3}{:<10}  {:<3}{:<3}{:<10}  {:<3}{:<3}{:<10}  {:<3}{:<3}{:<10}
{:<3}{:<3}{:<10}  {:<3}{:<3}{:<10}  {:<3}{:<3}{:<10}  {:<3}{:<3}{:<10}
""".format('A','α', 'alpha',    'B','β','beta',   'Γ','γ','gamma',      'Δ','δ', 'delta',
			'Ε','ε','epsilon',  'Z','ζ','zeta',   'H','η','eta',        'Θ','θ','theta',
			'I','ι','iota',     'K','κ','kappa',  '∧','λ','lambda',     'M','μ','mu',
			'N','ν','nu',       'Ξ','ξ','xi',     'Ο','ο','omicron',    '∏','π','pi',
			'P','ρ','rho',      'Σ','σ','sigma',  'T','τ','tau',        'Y','υ','upsilon',
			'Φ','φ','phi',      'X','x','chi',    'Ψ','ψ','psi',        'Ω','ω','omega')


order = input(" -> ").strip()
while order != 'exit':
	if order == "greek":
		print(greek)

	elif order == 'amino':
		print(amino_infomation)

	else:
		para = ' ' 
		if ' ' in order:
			order_list = order.split()
			order = order_list[0].strip()
			para = order_list[1].strip().strip('-')

		if order in ' '.join(amino_data):
			for item in amino_data:
				if order in [i.strip() for i in item.split(';')]:
					print(item)
					if para == 'p':
						im = Image.open('C:/Users/hhhhh/Desktop/databank/Enviroment/amino/' + item[0] + '.png')
						im.show()

		else:
			print('no results, check it! ')

	order = input(" -> ").strip()


