from flask import Flask, render_template, redirect, url_for, request, flash
from flask_bootstrap import Bootstrap
from flask_sqlalchemy import SQLAlchemy
from sqlalchemy.exc import SQLAlchemyError
from flask_wtf import FlaskForm
from wtforms import StringField, SubmitField, SelectField
from wtforms.validators import DataRequired
import requests
import re
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import plotly
import plotly.express as px
import json


app = Flask(__name__)
app.config['SECRET_KEY'] = '8BYkEfBA6O6donzWlSihBXox7C0sKR6b'
Bootstrap(app)
app.config['SQLALCHEMY_DATABASE_URI'] = "sqlite:///genes.db"
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
db = SQLAlchemy(app)


uniprot_search_url = "https://rest.uniprot.org/uniprotkb/search"
uniprot_fasta_url = "https://rest.uniprot.org/uniprotkb/"


class SearchForm(FlaskForm):
    search = StringField('Name of organism or protein:', validators=[DataRequired()])
    submit = SubmitField('Search')


class SelectForm(FlaskForm):
    first_gene = SelectField('First protein:', validators=[DataRequired()])
    second_gene = SelectField('Second protein:', validators=[DataRequired()])
    submit = SubmitField('Compare')


class Gene(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    uniprot_id = db.Column(db.String(250), unique=True, nullable=False)
    gene_name = db.Column(db.String(250), nullable=False)
    organism_name = db.Column(db.String(250), nullable=False)
    protein_name = db.Column(db.String(250), nullable=False)
    fasta = db.Column(db.String(3000), nullable=False)


# with app.app_context():
#     db.create_all()

def get_sequence(gene_id):
    requested_gene = Gene.query.get(gene_id)
    requested_sequence = requested_gene.fasta
    requested_sequence_cleared = re.split(r'SV=.\n', requested_sequence)[1].replace("\n", "").replace(" ", "")
    all_aminoacids = len(requested_sequence_cleared)
    protein_analysis = ProteinAnalysis(requested_sequence_cleared)
    aminoacid_count = protein_analysis.count_amino_acids()
    keys = list(aminoacid_count.keys())
    values = list(aminoacid_count.values())
    data_df = {
        'Aminoacid': keys,
        'Count': values
    }
    return [requested_gene, all_aminoacids, data_df]


@app.route('/')
def index():
    return render_template("index.html")


@app.route('/search', methods=["GET", "POST"])
def search():
    form = SearchForm()
    if form.validate_on_submit():
        query = form.search.data
        response = requests.get(uniprot_search_url, params={"query": query,
                                                            "format": "json",
                                                            "fields": "accession,id,organism_name,protein_name",
                                                            "size": 15
                                                            })
        data = response.json()
        amount_of_results = range(len(data["results"]))
        return render_template("select.html", results=data, amount=amount_of_results)
    return render_template("search.html", form=form)


@app.route("/atd")
def add_to_database():
    uniprot_id = request.args.get('uniprot_id')
    gene_name = request.args.get('gene_name')
    organism_name = request.args.get('organism_name')
    protein_name = request.args.get('protein_name')
    fasta_sequence = requests.get(f"{uniprot_fasta_url}{uniprot_id}.fasta").text
    new_gene = Gene(
        uniprot_id=uniprot_id,
        gene_name=gene_name,
        organism_name=organism_name,
        protein_name=protein_name,
        fasta=fasta_sequence,
    )
    try:
        db.session.add(new_gene)
        db.session.commit()
    except SQLAlchemyError as e:
        error = str(e.__dict__['orig'])
        flash(f"Something went wrong: {error}. Please try again.", "alert alert-danger")
        return redirect(url_for('search'))
    else:
        db.session.add(new_gene)
        db.session.commit()
        flash("Successfully added to database!", "alert alert-success")
        return redirect(url_for('search'))


@app.route("/database")
def database():
    all_genes = db.session.query(Gene).order_by('gene_name').all()
    return render_template('database.html', all_genes=all_genes)


@app.route("/dashboard/<int:gene_id>")
def dashboard(gene_id):
    sequence = get_sequence(gene_id)
    df = pd.DataFrame(sequence[2])
    fig = px.bar(df, x='Aminoacid', y='Count', title=f'Chart representing number of individual aminoacids in {sequence[0].gene_name} '
                                                     f'sequence, consisting of {sequence[1]} aminoacids')
    fig_json = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    return render_template('dashboard.html', gene=sequence[0], fig=fig_json)


@app.route("/compare", methods=['GET', 'POST'])
def compare():
    all_genes = db.session.query(Gene).order_by('protein_name').all()
    protein_names = [gene.protein_name for gene in all_genes]
    form = SelectForm()
    form.first_gene.choices = protein_names
    form.second_gene.choices = protein_names
    if form.validate_on_submit():
        first_gene = form.first_gene.data
        second_gene = form.second_gene.data
        first_gene_id = Gene.query.filter_by(protein_name=first_gene).first().id
        second_gene_id = Gene.query.filter_by(protein_name=second_gene).first().id
        return redirect(url_for('compare_dashboard', first_gene_id=first_gene_id, second_gene_id=second_gene_id))
    return render_template('compare.html', form=form)


@app.route("/compare-dashboard")
def compare_dashboard():
    first_gene_id = request.args.get('first_gene_id')
    second_gene_id = request.args.get('second_gene_id')
    first_sequence = get_sequence(first_gene_id)
    second_sequence = get_sequence(second_gene_id)
    first_gene_name = f"{first_sequence[0].gene_name} {first_sequence[0].protein_name}"
    second_gene_name = f"{second_sequence[0].gene_name} {second_sequence[0].protein_name}"
    first_gene_dict = first_sequence[2]
    second_gene_dict = second_sequence[2]
    data_df = {
        'Aminoacid': first_gene_dict['Aminoacid'] + second_gene_dict['Aminoacid'],
        'Count': first_gene_dict['Count'] + second_gene_dict['Count'],
        'Gene': [first_gene_name for a in first_gene_dict['Aminoacid']] +
                [second_gene_name for a in second_gene_dict['Aminoacid']]
    }
    df = pd.DataFrame(data_df)
    fig = px.bar(df, x='Aminoacid', y='Count', color='Gene', barmode='group',
                 title=f'Numbers of individual aminoacids in {first_gene_name} and {second_gene_name}')
    fig_json = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    return render_template('compare-dashboard.html', first_gene=first_sequence[0], second_gene=second_sequence[0],
                           fig=fig_json, first_gene_length=first_sequence[1], second_gene_length=second_sequence[1])

@app.route("/delete", methods=['GET', 'POST'])
def delete():
    gene_id = request.args.get('gene_id')
    gene = Gene.query.get(gene_id)
    try:
        db.session.delete(gene)
        db.session.commit()
    except SQLAlchemyError as e:
        error = str(e.__dict__['orig'])
        flash(f"Something went wrong: {error}. Please try again.", "alert alert-danger")
        return redirect(url_for('search'))
    else:
        db.session.delete(gene)
        db.session.commit()
        flash("Successfully deleted from database!", "alert alert-success")
        return redirect(url_for('database'))
    return redirect(url_for('database'))


if __name__ == '__main__':
    app.run(debug=True)
