#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sqlalchemy

print("The version of SQLAlchemy is " + str(sqlalchemy.__version__) )
print("-----------------")

from sqlalchemy import create_engine
engine = create_engine('sqlite:///:memory:', echo=True)



from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()


from sqlalchemy import Column, Integer, String
class User(Base):
     __tablename__ = 'users'
     id = Column(Integer, primary_key=True) # <-- define the primary key
     name = Column(String) # Users familiar with the syntax of CREATE TABLE may notice that the VARCHAR columns were generated without a length; on SQLite and Postgresql, this
     # is a valid datatype, but on others, it’s not allowed. So if running this tutorial on one of those databases, and you wish to use SQLAlchemy to issue CREATE TABLE, a “length” may be provided to the String type as below: Column(String(50))
     fullname = Column(String)
     password = Column(String)

     def __init__(self, name, fullname, password):
         self.name = name
         self.fullname = fullname
         self.password = password

     def __repr__(self):
        return "<User('%s','%s', '%s')>" % (self.name, self.fullname, self.password)



Base.metadata.create_all(engine)


ed_user = User('ed', 'Ed Jones', 'edspassword')



class User2(Base):
    __tablename__ = 'users2'
    id = Column(Integer, primary_key=True)
    name = Column(String)
    fullname = Column(String)
    password = Column(String)

u2 = User2(name='ed', fullname='Ed Jones', password='foobar')



from sqlalchemy.orm import sessionmaker
Session = sessionmaker(bind=engine)








#Session = sessionmaker()



session = Session()

session.add(ed_user)

session.add_all([
     User('wendy', 'Wendy Williams', 'foobar'),
     User('mary', 'Mary Contrary', 'xxg527'),
     User('fred', 'Fred Flinstone', 'blah')])


session.dirty
session.new

session.commit()
#session.add(u2)
