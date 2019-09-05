import React from 'react';

export default class Form extends React.Component {
    constructor() {
        super();
        this.handleSubmit = this.handleSubmit.bind(this);

    }

    handleSubmit(event) {
        event.preventDefault();
        const data = new FormData(event.target);
        fetch('/react-scatter-board', {
            method: 'POST',
            body: data
        })
    }

    render () {
        return (
            <div>
                <form onSubmit={this.handleSubmit}>
                    <input type="text" name="drug" id="drug"/>
                    <button>Submit</button>
                </form>
                
           </div>
        )
    }
}