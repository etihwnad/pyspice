
from pyparsing import *

title = StringStart() + SkipTo(LineEnd())
comment = StringStart() + Literal('*') + SkipTo(LineEnd())
element = StringStart() + SkipTo(LineEnd())
control = LineStart() + Literal('.') + SkipTo(LineEnd())
card = comment | element 
deck = title.setName('title') + OneOrMore(card)
deck = OneOrMore(card)

for line in open('hspiceFinal'):
    print card.parseString(line)
