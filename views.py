from flask import Flask, render_template, flash, redirect, session, url_for, request, g
from metmitdb import app
from sqlalchemy import or_
from sqlalchemy import func
from flask import render_template, Response, flash, redirect, url_for
from flask.ext.wtf import Form
from wtforms import StringField, BooleanField, TextAreaField
from wtforms.validators import DataRequired, Length
from models import *
from sequence_stats import get_sequence_stats
import collections
class SearchForm(Form):
    search = StringField('search', validators=[DataRequired()])

@app.before_request
def before_request():
    g.search_form = SearchForm()

@app.route('/')
@app.route('/index')
def index():
    return render_template('index.html')

@app.route('/search', methods=['POST'])
def search():
    if not g.search_form.validate_on_submit():
        return redirect(url_for('index'))
    return redirect(url_for('search_results', query=g.search_form.search.data))

@app.route('/search_results/<query>')
def search_results(query):
    results = db.session.query(Genome.name,Genome.ncbi_id).filter(or_(func.lower(Genome.name).contains(func.lower(query)),
                                            Genome.ncbi_id == query)).all()
    if (len(results)==0) or (len(results)>=20):
        flash('No genomes found for query %s. Try searching for something else' % query)
        return redirect(url_for('index'))
    return render_template('search_results.html',
                           query=query,
                           results=results)

@app.route('/single')
def single():
    return render_template('single_view.html')

@app.route('/ncbi_id/<ncbi_id>')
def genome_by_ncbi_id(ncbi_id):
    genome = Genome.query.filter_by(ncbi_id=ncbi_id).first()
    if genome == None:
        flash('mtDNA Genome with NCBI ID %s was not found. Try searching for something else' % ncbi_id)
        return redirect(url_for('index'))
    stat = collections.OrderedDict(sorted(get_sequence_stats(genome).items(),  reverse=True))
    return render_template('single_view.html',
                           genome=genome,
                           stat = stat)

@app.route('/get_fasta/<genome_id>.fasta')
def get_fasta(genome_id):
    genome = Genome.query.filter_by(id=genome_id).first()
    if genome.fasta == None:
        flash('Fasta for genome with genome id %s not found. Try searching for something else' % genome_id)
        return redirect(url_for('index'))
    else:
        return Response(genome.fasta, mimetype='text/plain')