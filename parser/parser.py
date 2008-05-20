
from pyparsing import *

title = Group(LineStart() + SkipTo(LineEnd()))
end = StringEnd()
#title = StringStart() + SkipTo(LineEnd())
emptyLine = LineStart() + SkipTo(LineEnd())
comment = Literal('*') + SkipTo(LineEnd())
elementStart = Word(alphanums) + SkipTo(LineEnd())
elementContinue = Literal('+') + SkipTo(LineEnd())
control = Literal('.') + SkipTo(LineEnd())

capacitor = (LineStart() + CaselessLiteral('c') + Word(alphanums) + Word(alphanums) + Word('alphanums') + Word(alphanums + '.') + ZeroOrMore(Word(alphanums + '.'))).setResultsName('capacitor')

element = capacitor ^ (elementStart + ZeroOrMore(elementContinue)).setResultsName('element')

card = Group(comment.setResultsName('comment') \
             | element \
             | emptyLine.setResultsName('empty') )

deck = title + ZeroOrMore(card) + end

infile = open('netlist').read()

result = deck.parseString(infile)

print result.dump()

for r in result:
    print r.dump()
#print result

#result = deck.parseFile(infile)



