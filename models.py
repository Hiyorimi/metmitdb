from metmitdb import db

class Genome(db.Model):
    __searchable__ = ['name','ncbi_id','accession_number']
    __tablename__ = 'genome'
    id = db.Column(db.Integer, primary_key=True)
    number = db.Column(db.Integer)
    ncbi_id = db.Column(db.Integer)
    link = db.Column(db.String)
    name = db.Column(db.String)
    bp_length = db.Column(db.Integer)
    accession_number = db.Column(db.String)
    gi = db.Column(db.Integer)
    fasta = db.Column(db.String)
    features = db.Column(db.Binary)

    def __repr__(self):
        return "<Genome(ncbi_id='%s')>" % (
        self.ncbi_id)
