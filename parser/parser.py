
from pyparsing import *
from decimal import Decimal

title = Group(LineStart() + SkipTo(LineEnd()))
end = StringEnd()
name = Word(alphas, alphanums+'_')
node = Word(alphanums, alphanums+'_')
#title = StringStart() + SkipTo(LineEnd())
emptyLine = LineStart() + SkipTo(LineEnd())
comment = Literal('*') + SkipTo(LineEnd())
elementStart = Word(alphanums) + SkipTo(LineEnd())
elementContinue = Literal('+').suppress() + SkipTo(LineEnd())
control = Literal('.') + SkipTo(LineEnd())

point = Literal('.')
e = CaselessLiteral('e')
fnumber = Combine( Word( "+-"+nums, nums ) + 
                   Optional( point + Optional( Word( nums ) ) ) +
                   Optional( e + Word( "+-"+nums, nums ) ) )
value = fnumber

parameter = Group(name.setResultsName('name')
                  + Literal('=').suppress()
                  + value.setResultsName('val').setParseAction(lambda s,l,t: Decimal(t[0]))
             )
            

capacitor = Group(name.setResultsName('name')
                  + node.setResultsName('n0')
                  + node.setResultsName('n1')
                  + value.setResultsName('value')
                  + Group(ZeroOrMore(parameter)).setResultsName('parm')
                 )

basicElement = (elementStart +
                ZeroOrMore(elementContinue)).setResultsName('basicElement')

element = capacitor ^ basicElement

card = (Group(comment.setResultsName('comment') \
             | element \
             | emptyLine.setResultsName('empty') )
       ).setResultsName('card')

deck = Group(title.setResultsName('title') + ZeroOrMore(card) + end)

infile = open('netlist').read()

result = deck.parseString(infile)

#print result.dump()

for r in result:
    print r.dump()
#print result

teststr = '''title line
C100 na nb 0.001 k1=1.0 k2=2.5
C203 z1 z2 5.6 x2 = 4.1'''

print
td = deck.parseString(teststr)
for r in td:
    pass
    #print r.dump()

print

i=1
for r in td:
    print i
    print r.title
    print r.dump()
    j=1
    for c in r.card:
        print i,j
        print c.dump()
        j += 1
    i += 1


